package Bio::Tools::Run::Phylo::ExaBayes;
use strict;
use version;
use Cwd;
use File::Spec;
use File::Temp 'tempfile';
use Bio::Phylo::Generator;
use Bio::Phylo::IO 'parse';
use Bio::Tools::Run::Phylo::PhyloBase;
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::Util::Logger;
use Bio::Phylo::PhyLoTA::Service::TreeService;

use base qw(Bio::Tools::Run::Phylo::PhyloBase);

my $treeservice = Bio::Phylo::PhyLoTA::Service::TreeService->new;
my $log = Bio::Phylo::Util::Logger->new;

# program name required by PhyloBase
our $PROGRAM_NAME = 'exabayes'; 

our @ExaBayes_PARAMS = (
        'f',        # alnFile:    alignment file (binary and created by parser or plain-text phylip)
        's',        # seed:       master seed for the MCMC
        'n',        # ruid:       a run id
        'r',        # id:         restart from checkpoint id
        'q',        # modelfile:  RAxML-style model file
        't',        # treeFile:   file containing starting trees (newick format) for chains
        'm',        # model:      data type for single partition non-binary alignment file (DNA | PROT)
        'c',        # confFile:   file configuring your ExaBayes run
        'w',        # dir:        working directory for output files
        'R',        # num:        number of runs (i.e., independent chains) to be executed in parallel
        'C',        # num:        number of chains (i.e., coupled chains) to be executed in parallel
        'M'         # mode:       memory versus runtime trade-off, value 0 (fastest) to 3 (most memory efficient)        
    );

our @ExaBayes_CONFIGFILE_PARAMS = (
        # RUN parameters
        'numRuns',           # number of independant runs
        'numGen',            # number of generations to run. Defaults to 1e6
        'parsimonyStart',    # if 1, start from most parsimoneous tree
        #MCMC parameters
        'numCoupledChains'   # number of chains per independent run
    );

our @Consense_PARAMS = (
        # The following command line parameters are not for "ExaBayes" per se 
        # but for the tool "Consense" which is also handled in this class.
        't',        # thresh:     threshold for the consenus tree 
                    #             values between 50 (majority rule) and 100 (strict) or MRE (the greedily refined MR consensus).  
                    #             Default: MRE
        'b',        # relBurnin:  proportion of trees to discard as burn-in (from start). Default: 0.25        
    );

our @Misc_PARAMS = (
        # The following parameters are artificial, not a 'real' parameter of ExaBayes nor Consense
        'outfile_name',        # file name of the inferred consensus trees
        'outfile_format',       # format of output file. Either 'nexus' or 'newick' 
        'consense_bin',         # location of consense binary 
    );


our @ExaBayes_SWITCHES = (
        'v',        # version
        'z',        # quiet mode
        'd',        # dry run
        'Q',        # per-partition data distribution
        'S'         # try to save memory using the SEV-technique for gap columns       
    );


# Create function aliases for (one letter) command line arguments and switches

=item work_dir

Getter/setter. Intermediate files will be written here.

=cut

*work_dir = *w;

=item run_id

Getter/setter for an ID string for this run. by default based on the PID. this only has
to be unique during the runtime of this process, all files that have the run ID in them
are cleaned up afterwards.

=cut

*run_id = *n;

=item seed

Getter/setter for master seed for MCMC.

=cut

*seed = *s;

=item burnin_fraction

Getter/setter for fraction of burn in (used in the consense program).

=cut

*burnin_fraction = *b;

=item config_file

Getter/setter for a config file for the ExaBayes run.

=cut

*config_file = *c;

=item model

Getter/setter for Exabayes model (DNA | Protein)

=cut

*model = *m;



sub new {
        my ( $class, @args ) = @_;
        my $self = $class->SUPER::new(@args);
        $self->_set_from_args(
                \@args,
                '-methods' => [ @ExaBayes_PARAMS, @ExaBayes_SWITCHES, @ExaBayes_CONFIGFILE_PARAMS, @Misc_PARAMS, @Consense_PARAMS ],
                '-create'  => 1,
            );
        my ($out) = $self->SUPER::_rearrange( [qw(OUTFILE_NAME)], @args );
        $self->outfile_name( $out || '' );
        return $self;
}


sub run {
	my $self = shift;
        my $ret;
        my %args = @_;
        my $phylip  = $args{'-phylip'} || die "Need -phylip arg";
        my $intree  = $args{'-intree'} || die "Need -intree arg";
        
        my $binary = $self->run_id . '-dat' ;
        $binary = $treeservice->make_phylip_binary( $phylip, $binary, $self->parser, $self->work_dir );
        $binary = File::Spec->catfile($self->work_dir, $binary);
        my @tipnames = $treeservice->read_tipnames( $phylip );
        
        #my $forest = Bio::Phylo::IO->parse(
        #        '-format' => 'newick',
        #        '-file'   => $intree,
        #    );
        #my $tree = $forest->first;
        #$tree -> keep_tips( [ @tipnames ] );
        #$tree -> resolve;
        #$tree -> remove_unbranched_internals;
        #$tree -> deroot;                                                                      
        #my $newtree = $intree."_new";
        #open my $fh, '>', $newtree or die $!;
        #close $fh;        
        #my $intree = $self->_make_intree($taxa, $tree, \@tipnames );
        
        # compose argument string: add MPI commands, if any
        my $string;
        if ( $self->mpirun && $self->nodes ) {
		$string = sprintf '%s -np %i ', $self->mpirun, $self->nodes;
	}
        ##TODO: fix issue with input tree
        my $newtree = 1;
        $string .= $self->executable . $self->_setparams($binary, $newtree);       
                
        $log->info("going to run '$string'");        
        system($string) and $self->warn("Couldn't run ExaBayes: $?");
        
        # build command string for running 'consense'
        $string = $self->consense_bin . $self->_set_consense_params;

        # add all topology files created by Exabayes
        $string .= " -f " . $self->work_dir . "/ExaBayes_topologies." . $self->run_id . ".*" . " -n " . $self->run_id;
        
        # run consense
        $log->info("Going to run $string");        
        my $output = system($string) and $self->warn("Couldn't run Consense: $?");
        
        # get format for output file, if not given, chose 'newick'
        my $format = $self->outfile_format | "newick";
        $format = ucfirst($format);
        # consense output file name is not controllable, therefore rename consense outfile
        my $filename = "ExaBayes_Consensus*".$format . "." . $self->run_id;
        
        if ( $self->outfile_name ) {
                my $mvstr = "mv " . $filename . " " . $self->outfile_name;
                system ($mvstr) and $self->warn("Couldn't run command $mvstr : $?");
                $ret = $self->outfile_name;
        } 
        else {
                $ret = $filename;
        }
        #return $self->_cleanup;
        return $ret;
}

=item mpirun

Getter/setter for the location of the mpirun program for parallelized runs. When unset,
no parallelization is attempted.

=cut

sub mpirun {
        my ( $self, $mpirun ) = @_;
        if ( $mpirun ) {
		$self->{'_mpirun'} = $mpirun;
	}
	return $self->{'_mpirun'};
}

=item nodes

Getter/setter for number of mpi nodes

=cut

sub nodes {
	my ( $self, $nodes ) = @_;
	if ( $nodes ) {
		$self->{'_nodes'} = $nodes;
	}
	return $self->{'_nodes'};
}

=item parser

Getter/setter for the location of the parser program (which creates a compressed, binary 
representation of the input alignment) that comes with ExaBayes. 

=cut

sub parser {
	my ( $self, $parser ) = @_;
	if ( $parser ) {
		$self->{'_parser'} = $parser;
	}
	return $self->{'_parser'} || 'parser';
}

sub _run_in_dir {
        my ( $self, $dir, $command ) = @_;
        #return value
        my $res;
        if ( $dir && $command ) {
                $log->info("Going to run command $command in $dir");
                my $curdir = getcwd;
                chdir $self->work_dir;	
                $res = system($command) and $self->warn("Couldn't run command $command : $?");                
                chdir $curdir;
        }
        else {
                $log->warn("Need 'directory' and 'command' arguments");
                $res = -1;
        }        
        return $res;
}

sub _write_config_file {
        my $self = shift;
        my $conffile = File::Spec->catfile( $self->work_dir, $self->run_id . '.nex' );
	open my $conffh, '>', $conffile or die $!;
        print $conffh "#NEXUS\n";
        print $conffh "begin run;\n";
        for my $attr (@ExaBayes_CONFIGFILE_PARAMS) {
                my $value = $self->$attr();
                next unless defined $value;
                print $conffh $attr . "\t" . $value . "\n";
        }
        print $conffh "end;\n";
        close $conffh;
        $self->config_file($conffile);
        return $self->config_file;
}

sub _make_intree {
	my ( $self, $tree, $tipnames ) = @_;
	my $treefile = File::Spec->catfile( $self->work_dir, $self->run_id . '.dnd' );
	open my $treefh, '>', $treefile or die $!;
	if ( $tree ) {
                $tree -> keep_tips( [@{$tipnames}] );
                $tree -> resolve;
                $tree -> remove_unbranched_internals;
                $tree -> deroot;                                                              
                print $treefh $tree->to_newick;
	}
	else {
	
		# no tree was given to the class. here we then simulate
		# a BS tree shape.
		my $gen = Bio::Phylo::Generator->new;
		$tree = $gen->gen_equiprobable( '-tips' => scalar @{$tipnames} )->first;
		my $i = 0;
		$tree->visit(sub{
			my $n = shift;
			if ( $n->is_terminal ) {
				$n->set_name( @{$tipnames}[$i++] );
			}
		});
                return $self->_make_intree( $tree, $tipnames );
	}
	return $treefile;
}



sub _setparams {
    my ( $self, $infile, $intree ) = @_;
    my $param_string = '';

    # add config file to parameters if not already existant
    if (! $self->config_file){
            $self->_write_config_file;
    }
    #$params_string .= ' -c ' . self->$config_file;

    # iterate over parameters and switches
    for my $attr (@ExaBayes_PARAMS) {
        my $value = $self->$attr();
        next unless defined $value;
        $param_string .= ' -' . $attr . ' ' . $value;
    }
    for my $attr (@ExaBayes_SWITCHES) {
        my $value = $self->$attr();
        next unless $value;
        $param_string .= ' -' . $attr;
    }
    
    # set file names to local
    my %path = ('-f' => $infile );  ###, '-t' => $intree  , '-n' => $self->outfile_name );
    while( my ( $param, $path ) = each %path ) {		
		$param_string .= " $param $path";
    }
    
    # hide stderr
    my $null = File::Spec->devnull;
    $param_string .= " > $null 2> $null" if $self->quiet() || $self->verbose < 0;

    return $param_string;
}

sub _set_consense_params {
        my $self = shift;
        my $param_string;
        for my $attr (@Consense_PARAMS) {
                        my $value = $self->$attr();
                        next unless defined $value;
                        $param_string .= ' -' . $attr . ' ' . $value;                        
        }
        # hide stderr
        my $null = File::Spec->devnull;
        $param_string .= " > $null 2> $null" if $self->quiet() || $self->verbose < 0;
        return $param_string;
}


##sub program_name { $PROGRAM_NAME }

=item program_dir

(no-op)

=cut

sub program_dir { undef }


=item version

Returns the version number of the ExaBayes executable.

=cut

sub version {
    my ($self) = @_;
    my $exe;
    return undef unless $exe = $self->executable;
    my $string = `$exe -v 2>&1`;
    $string =~ /ExaBayes, version (\d+\.\d+\.\d+)/;
    return version->parse($1) || undef;
}

1;
