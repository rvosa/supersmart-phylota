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

use Bio::Tools::Run::Phylo::ExaML;

use base qw(Bio::Tools::Run::Phylo::PhyloBase);

my $examl = Bio::Tools::Run::Phylo::ExaML->new;

our $PROGRAM_NAME = 'yggdrasil';
our $POSTPROCESS_PROGRAM_NAME = 'consense';

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
        'outfile_name',        #  file name of the inferred consensus trees
        'outfile_format'       #  format of output file. Either 'nexus' or 'newick' 
    );


our @ExaBayes_SWITCHES = (
        'v',        # version
        'z',        # quiet mode
        'd',        # dry run
        'Q',        # per-partition data distribution
        'S'         # try to save memory using the SEV-technique for gap columns       
    );


my $log = Bio::Phylo::Util::Logger->new;

##*work_dir = *w;

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    $self->_set_from_args(
        \@args,
        '-methods' => [ @ExaBayes_PARAMS, @ExaBayes_SWITCHES ],
        '-create'  => 1,
    );
    my ($out) = $self->SUPER::_rearrange( [qw(OUTFILE_NAME)], @args );
    $self->outfile_name( $out || '' );
    return $self;
}

#sub run {
#    my $self = shift;
#    my $phylip;
#    if (@_ == 1){
#	$phylip = shift;
#    } else {
#	die "Need phylip file";
#    }
#    print "PHYLIP : $phylip \n";

    # compile command string with executable and parameters
#    my $string = $self->executable . $self->_setparams($phylip, @ExaBayes_PARAMS, @ExaBayes_SWITCHES);       
#    print "String : $string \n";
#    
#    my $curdir = getcwd;
#    $log->info("going to run '$string' inside ".$self->w);
#    print "going to run '$string' inside ".$self->w."\n";
#    system($string) and $self->warn("Couldn't run ExaBayes: $?");

    # compile command string for building the conensus tree
    ##my $post_string = $POSTPROCESS_PROGRAM_NAME . 
    
#}

=item 
sub prepare_tree {
        my $self = shift;
        my ($intree, $phylip) = @_;
        
        my @tipnames = read_tipnames( $phylip );
        my $tree = parse(
                '-format' => 'newick',
                '-file'   => $intree,
            );
        $tree -> keep_tips( @tipnames );
        $tree -> resolve;
        $tree -> remove_unbranched_internals;
        $tree -> deroot;                                      
        return($examl->_make_intree( $taxa, $tree ));
} 
=cut

=item c
sub run {
	my $self = shift;	
        my %args = @_;
        my $phylip  = $args{'-phylip'} || die "Need -phylip arg";
        my $intree  = $args{'-intree'} || die "Need -intree arg";
        my $tree = parse(
                '-format' => 'newick',
                '-file'   => $intree,
            );                
        my $binary = $examl->_make_binary( $phylip );               
        my $string = $self->executable . $self->_setparams($binary, $intree);       
        print "String : $string \n";
        
        my $curdir = getcwd;
        $log->info("going to run '$string' inside ".$self->w);
        print "going to run '$string' inside ".$self->w."\n";
        system($string) and $self->warn("Couldn't run ExaBayes: $?");
} 
=cut 

sub run {
	my $self = shift;
 	my ($project, $phylip);# $phylip, $intree );
	
	# one argument means it's a nexml file
	if ( @_ == 1 ) {
                print "HIER\n";
                my $nexml = shift;	
		# read nexml input
		$project = parse(
			'-format'     => 'nexml',
			'-file'       => $nexml,
			'-as_project' => 1,
		);	
	}
	else {
		my %args = @_;
                $project = parse(
			'-format'     => 'newick',
			'-file'       => $args{'-intree'},
			'-as_project' => 1,
		);
                $phylip = $args{'-phylip'}
	}
        print "Project ".ref($project)."\n";
        my @matrix = @{ $project->get_items(_MATRIX_) };
        my ($taxa) = @{ $project->get_items(_TAXA_) };
        print "Matrix ".ref(@matrix[0])."\n";        
        print "Taxa ".ref($taxa)."\n";        
        
        if (! $phylip){                                
                $phylip = $examl->_make_phylip( $taxa, @matrix );
        }
        print "Phylip : $phylip \n";
        my $binary = $examl->_make_binary( $phylip );
        my @tipnames = read_tipnames( $phylip );

        
        #for my $ti (@tipnames){
        #        print "Ti : $ti \n";                
        #}
        my ($tree) = @{ $project->get_items(_TREE_) }; 
        my $intree = $self->_make_intree($taxa, $tree, \@tipnames );
        print "Tree ".ref($tree)."\n";

        my $string = $self->executable . $self->_setparams($binary, $intree);       
        print "String : $string \n";
        
        my $curdir = getcwd;
        $log->info("going to run '$string' inside ".$self->w);
        print "going to run '$string' inside ".$self->w."\n";
        system($string) and $self->warn("Couldn't run ExaBayes: $?");
}


=item parser

Getter/setter for the location of the parser program (which creates a compressed, binary 
representation of the input alignment) that comes with examl. By default this is found
on the PATH.

=cut

sub parser {
	my ( $self, $parser ) = @_;
	if ( $parser ) {
		$self->{'_parser'} = $parser;
	}
	return $self->{'_parser'} || 'parser';
}


sub _make_intree {
	my ( $self, $taxa, $tree, $tipnames ) = @_;
	my $treefile = File::Spec->catfile( $self->w, $self->n . '.dnd' );
	open my $treefh, '>', $treefile or die $!;
	if ( $tree ) {
                ##for my $ti (@{$tipnames}){
                        ##print "Ti2:$ti: \n";
                ##}                                        
                $tree -> keep_tips( [@{$tipnames}] );
                $tree -> resolve;
                $tree -> remove_unbranched_internals;
                $tree -> deroot;                                                              
                print $treefh $tree->to_newick;
	}
	else {
	
		# no tree was given in the nexml file. here we then simulate
		# a BS tree shape.
		my $gen = Bio::Phylo::Generator->new;
		$tree = $gen->gen_equiprobable( '-tips' => scalar @{$tipnames} )->first;
		my $i = 0;
		$tree->visit(sub{
			my $n = shift;
			if ( $n->is_terminal ) {
				$n->set_name( @{$tipnames}[$i++]);
			}
		});
                return $self->_make_intree( $taxa, $tree, $tipnames);
	}
	return $treefile;
}


# reads the supermatrix (phylip format) and returns the tip names from it
sub read_tipnames {
	my $supermatrix = shift;
	my $ntax = 0;
	my $line = 0;
	my @result;	
        print "Supermatrix : $supermatrix \n";
        open my $fh, '<', $supermatrix or die $!;
	LINE: while(<$fh>) {
		chomp;
		my $word;
		if ( /^(\S+)/ ) {
			$word = $1;
		}
		if ( not $ntax ) {
			$ntax = $word;
			##$logger->debug("$supermatrix has $ntax taxa");
			next LINE;
		}
		push @result, $word;
		##$logger->debug("adding taxon $word");
		last LINE if ++$line == $ntax;
	}
	return @result;
}



sub _setparams {
    my ( $self, $infile, $intree ) = @_;
    my $param_string = '';

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

    # Set default output file if no explicit output file has been given
#    if ( ! $self->outfile_name ) {
#        my ( $tfh, $outfile ) = $self->io->tempfile( '-dir' => $self->work_dir );
#        close $tfh;
#        undef $tfh;
#        $self->outfile_name($outfile);
#    }
    
    # set file names to local
    my %path = ('-f' => $infile,  '-t' => $intree ); ##, '-n' => $self->outfile_name );
	while( my ( $param, $path ) = each %path ) {		
		$param_string .= " $param $path";
	}
    
    # hide stderr
    my $null = File::Spec->devnull;
    $param_string .= " > $null 2> $null" if $self->quiet() || $self->verbose < 0;

    return $param_string;
}



#sub _setparams {
#    my ( $self, $infile, @params, @switches, $binary, $intree) = @_;
#    my $param_string = '';

    # iterate over parameters and switches
#    for my $attr (@params) {
#        my $value = $self->$attr();
#        next unless defined $value;
#        $param_string .= ' -' . $attr . ' ' . $value;
#    }
#    for my $attr (@switches) {
#        my $value = $self->$attr();
#        next unless $value;
#        $param_string .= ' -' . $attr;
#    }

    # Set default output file if no explicit output file has been given
#    if ( ! $self->outfile_name ) {
#            my ( $tfh, $outfile ) = $self->io->tempfile( '-dir' => $self->work_dir );
#            close $tfh;
#            undef $tfh;
#            $self->outfile_name($outfile);
#    }
    
    # set file names to local
    #my %path = ( '-t' => $intree, '-s' => $infile, '-n' => $self->outfile_name );
#    my %path = ('-f' => $infile);
#    while( my ( $param, $path ) = each %path ) {
#            $param_string .= " $param $path";
#    }
    
    # hide stderr
#    my $null = File::Spec->devnull;
#    $param_string .= " > $null 2> $null" if $self->quiet() || $self->verbose < 0;

#    return $param_string;
#}


sub program_name { $PROGRAM_NAME }

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
