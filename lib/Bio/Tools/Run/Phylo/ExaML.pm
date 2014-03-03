package Bio::Tools::Run::Phylo::ExaML;
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
use base qw(Bio::Tools::Run::Phylo::PhyloBase);

=head1 NAME

Bio::Tools::Run::Phylo::ExaML - tree inference wrapper for ExaML

=head1 SYNOPSYS

 use Bio::Tools::Run::Phylo::ExaML;
 my $examl Bio::Tools::Run::Phylo::ExaML->new;
 $examl->outfile_name('out.tre');
 $examl->work_dir('.');
 $examl->m('GAMMA');
 $examl->quiet(1);
 $examl->run('infile.phy');

=head1 DESCRIPTION

Infers a large phylogeny using the ExaML program. Input files must be in PHYLIP format,
(so far this has only been tested with interleaved PHYLIP), output file contains an
unrooted tree in Newick format.

=head1 METHODS

B<Note:> in addition to the methods described below, all the arguments described for the
constructor can be used as methods as well. For example, the argument C<-B> becomes:

 $examl->B($num);

=over

=cut

our $PROGRAM_NAME = 'examl';
our @ExaML_PARAMS = qw(B c e f i m);
our @ExaML_SWITCHES = qw(a D M Q S);
my $log = Bio::Phylo::Util::Logger->new;

=item new

The constructor takes the following optional, named arguments that require appropriate 
values:

 '-B' # specify the number of best ML trees to save and print to file
 '-c' # Specify number of distinct rate catgories for ExaML when modelOfEvolution
        is set to GTRPSR. Individual per-site rates are categorized into 
        numberOfCategories rate categories to accelerate computations. DEFAULT: 25
 '-e' # set model optimization precision in log likelihood units for final
        optimization of model parameters. DEFAULT: 0.1
 '-f' # select algorithm: 
 			'-f' => 'd': new rapid hill-climbing DEFAULT: ON
            '-f' => 'o': old and slower rapid hill-climbing without heuristic cutoff
 '-i' # Initial rearrangement setting for the subsequent application of topological 
        changes phase
 '-m' # Model of rate heterogeneity:
			'-m' => 'PSR' for the per-site rate category model (this used to be called 
			        CAT in RAxML)
			'-m' => 'GAMMA' for the gamma model of rate heterogeneity with 4 discrete 
         	        rates

In addition, the following optional, named arguments that require booleans are available:

 '-a' # use the median for the discrete approximation of the GAMMA model of rate 
        heterogeneity. DEFAULT: OFF
 '-D' # ML search convergence criterion. This will break off ML searches if the relative 
        Robinson-Foulds distance between the trees obtained from two consecutive lazy SPR 
        cycles is smaller or equal to 1%. Usage recommended for very large datasets in 
        terms of taxa. On trees with more than 500 taxa this will yield execution time 
        improvements of approximately 50% while yielding only slightly worse trees.
        DEFAULT: OFF
 '-M' # Switch on estimation of individual per-partition branch lengths. Only has effect 
        when used in combination with "-q". Branch lengths for individual partitions will 
        be printed to separate files. A weighted average of the branch lengths is computed 
        by using the respective partition lengths.
 '-Q' # Enable alternative data/load distribution algorithm for datasets with many 
        partitions. In particular under PSR this can lead to parallel performance 
        improvements of up to factor 10!
 '-S' # turn on memory saving option for gappy multi-gene alignments. For large and gappy 
        datasets specify -S to save memory. This will produce slightly different 
        likelihood values, may be a bit slower but can reduce memory consumption from 
        70GB to 19GB on very large and gappy datasets.

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    $self->_set_from_args(
        \@args,
        '-methods' => [ @ExaML_PARAMS, @ExaML_SWITCHES ],
        '-create'  => 1,
    );
    my ($out) = $self->SUPER::_rearrange( [qw(OUTFILE_NAME)], @args );
    $self->outfile_name( $out || '' );
    return $self;
}

=item program_name

Returns the local name of the wrapped program, i.e. 'examl'.

=cut

sub program_name { $PROGRAM_NAME }

=item program_dir

(no-op)

=cut

sub program_dir { undef }

=item run_id

Getter/setter for an ID string for this run. by default based on the PID. this only has
to be unique during the runtime of this process, all files that have the run ID in them
are cleaned up afterwards.

=cut

sub run_id {
	my $self = shift;
	if ( @_ ) {
		$self->{'_run_id'} = shift;
	}
	return $self->{'_run_id'} || "examl-run-$$";
}

=item work_dir

Getter/setter. Intermediate files will be written here.

=cut

sub work_dir {
	my $self = shift;
	if ( @_ ) {
		$self->{'_workdir'} = shift;
	}
	return $self->{'_workdir'} || '.';
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

=item version

Returns the version number of the ExaML executable.

=cut

sub version {
    my ($self) = @_;
    my $exe;
    return undef unless $exe = $self->executable;
    my $string = `$exe -v 2>&1`;
    $string =~ /ExaML version (\d+\.\d+\.\d+)/;
    return version->parse($1) || undef;
}

=item run

Runs the analysis. Given a single argument, this is assumed to be the location of a
NeXML file that contains both one or more alignments and a starting tree. Given the
named arguments C<-phylip> and C<-intree> (both are required in combination), uses
these as the data and starting tree.

=back

=cut

sub run {
	my $self = shift;
	my ( $phylip, $intree );
	
	# one argument means it's a nexml file
	if ( @_ == 1 ) {
		my $nexml = shift;
	
		# read nexml input
		my $project = parse(
			'-format'     => 'nexml',
			'-file'       => $nexml,
			'-as_project' => 1,
		);
		my ($tree) = @{ $project->get_items(_TREE_) };
		my @matrix = @{ $project->get_items(_MATRIX_) };
		my ($taxa) = @{ $project->get_items(_TAXA_) };
	
		# create input files
		$phylip = $self->_make_phylip( $taxa, @matrix );
		$intree = $self->_make_intree( $taxa, $tree );		
	}
	else {
		my %args = @_;
		$phylip  = $args{'-phylip'} || die "Need -phylip arg";
		$intree  = $args{'-intree'} || die "Need -intree arg";
	}
	my $binary = $self->_make_binary( $phylip );
	
	# compose argument string: add MPI commands, if any
	my $string;
	if ( $self->mpirun && $self->nodes ) {
		$string = sprintf '%s -np %i ', $self->mpirun, $self->nodes;
	}
	
	# add executable and parameters
	$string .= $self->executable . $self->_setparams($binary,$intree);
	
	# examl wants to run inside the dir with data
	my $curdir = getcwd;
	chdir $self->work_dir;	
	$log->info("going to run '$string' inside ".$self->work_dir);
	system($string) and $self->warn("Couldn't run ExaML: $?");
	chdir $curdir;
	
	# remove cruft
	return $self->_cleanup;
}

sub _cleanup {
	my $self = shift;
	my $dir  = $self->work_dir;
	my ( $outv, $outd, $out ) = File::Spec->splitpath($self->outfile_name);
	my $run  = $self->run_id;
	opendir my $dh, $dir or die $!;
	while( my $entry = readdir $dh ) {
		if ( $entry =~ /^ExaML_binaryCheckpoint\.${out}_\d+$/ ) {
			unlink "${dir}/${entry}";
		}
		elsif ( $entry =~ /^ExaML_(?:info|log).${out}$/ ) {
			unlink "${dir}/${entry}";		
		}
		elsif ( $entry =~ /^${run}-dat\.binary$/ ) {
			unlink "${dir}/${entry}";		
		}
		elsif ( $entry =~ /^${run}\.(?:dnd|phy)$/ ) {
			unlink "${dir}/${entry}";		
		}
		elsif ( $entry =~ /^RAxML_info\.${run}-dat$/ ) {
			unlink "${dir}/${entry}";		
		}
	}
	rename "${dir}/ExaML_result\.${out}", $self->outfile_name;
	return $self->outfile_name;
}

sub _make_binary {
	my ( $self, $phylip ) = @_;
	my $binfile = File::Spec->catfile( $self->work_dir, $self->run_id . '-dat' );
	$log->info("going to make binary representation of $phylip => $binfile");	
	my ( $binvolume, $bindirectories, $binbase ) = File::Spec->splitpath( $binfile );
	my ( $phylipvolume, $phylipdirectories, $phylipbase ) = File::Spec->splitpath( $phylip );	
	my $curdir = getcwd;
	chdir $self->work_dir;
	my @command = ( $self->parser, 
		'-m' => 'DNA', 
		'-s' => $phylipbase, 
		'-n' => $binbase,
		'>'  => File::Spec->devnull,		
		'2>' => File::Spec->devnull,
	);
	my $string = join ' ', @command;
	$log->info("going to run '$string' inside ".$self->work_dir);
	system($string) and $self->warn("Couldn't create $binfile: $?");
	chdir $curdir;
	return "${binfile}.binary";
}

sub _make_intree {
	my ( $self, $taxa, $tree ) = @_;
	my $treefile = File::Spec->catfile( $self->work_dir, $self->run_id . '.dnd' );
	open my $treefh, '>', $treefile or die $!;
	if ( $tree ) {
		print $treefh $tree->to_newick;
	}
	else {
	
		# no tree was given in the nexml file. here we then simulate
		# a BS tree shape.
		my $gen = Bio::Phylo::Generator->new;
		$tree = $gen->gen_equiprobable( '-tips' => $taxa->get_ntax )->first;
		my $i = 0;
		$tree->visit(sub{
			my $n = shift;
			if ( $n->is_terminal ) {
				$n->set_name( $taxa->get_by_index($i++)->get_name );
			}
		});
		return $self->_make_intree( $taxa, $tree );
	}
	return $treefile;
}

sub _make_phylip {
	my ( $self, $taxa, @matrix ) = @_;
	
	# create phylip file for parser
	my $phylipfile = File::Spec->catfile( $self->work_dir, $self->run_id . '.phy' );
	open my $phylipfh, '>', $phylipfile or die $!;
	my %nchar_for_matrix;
	my $ntax  = $taxa->get_ntax;
	my $nchar = 0;
	for my $m ( @matrix ) {
		my $mid = $m->get_id;
		$nchar_for_matrix{$mid} = $m->get_nchar;
		$nchar += $nchar_for_matrix{$mid};
	}
	print $phylipfh $ntax, ' ', $nchar, "\n";
	$taxa->visit(sub{
		my $t = shift;
		my @d = @{ $t->get_data };
		print $phylipfh $t->get_name, ' ';
		for my $m ( @matrix ) {
			my $mid = $m->get_id;
			my ($row) = grep { $_->get_matrix->get_id == $mid } @d;
			if ( $row ) {
				print $phylipfh $row->get_char;
			}
			else {
				print $phylipfh '?' x $nchar_for_matrix{$mid};
			}
		}
		print $phylipfh "\n";
	});	
	return $phylipfile;
}

sub _setparams {
    my ( $self, $infile, $intree ) = @_;
    my $param_string = '';

	# iterate over parameters and switches
    for my $attr (@ExaML_PARAMS) {
        my $value = $self->$attr();
        next unless defined $value;
        $param_string .= ' -' . $attr . ' ' . $value;
    }
    for my $attr (@ExaML_SWITCHES) {
        my $value = $self->$attr();
        next unless $value;
        $param_string .= ' -' . $attr;
    }

    # Set default output file if no explicit output file has been given
    if ( ! $self->outfile_name ) {
        my ( $tfh, $outfile ) = $self->io->tempfile( '-dir' => $self->work_dir );
        close $tfh;
        undef $tfh;
        $self->outfile_name($outfile);
    }
    
    # set file names to local
    my %path = ( '-t' => $intree, '-s' => $infile, '-n' => $self->outfile_name );
	while( my ( $param, $path ) = each %path ) {
		my ( $volume, $directories, $file ) = File::Spec->splitpath( $path );
		$param_string .= " $param $file";
	}
    
    # hide stderr
    my $null = File::Spec->devnull;
    $param_string .= " > $null 2> $null" if $self->quiet() || $self->verbose < 0;

    return $param_string;
}

1;