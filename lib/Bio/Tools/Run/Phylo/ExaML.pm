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
use base qw(Bio::Tools::Run::Phylo::PhyloBase);

our $PROGRAM_NAME = 'examl';
our @ExaML_PARAMS = qw(B c e f i m);
our @ExaML_SWITCHES = qw(a D M Q S);

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

sub program_name { $PROGRAM_NAME }

sub program_dir { undef }

sub run_id {
	my $self = shift;
	if ( @_ ) {
		$self->{'_run_id'} = shift;
	}
	return $self->{'_run_id'} || "examl-run-$$";
}

# intermediate files will be written here
sub work_dir {
	my $self = shift;
	if ( @_ ) {
		$self->{'_workdir'} = shift;
	}
	return $self->{'_workdir'} || '.';
}

# getter/setter for the parser program (which creates a compressed, binary representation
# of the input alignment) that comes with examl
sub parser {
	my ( $self, $parser ) = @_;
	if ( $parser ) {
		$self->{'_parser'} = $parser;
	}
	return $self->{'_parser'} || 'parser';
}

sub version {
    my ($self) = @_;
    my $exe;
    return undef unless $exe = $self->executable;
    my $string = `$exe -v 2>&1`;
    $string =~ /ExaML version (\d+\.\d+\.\d+)/;
    return version->parse($1) || undef;
}

sub run {
	my ( $self, $nexml ) = @_;
	
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
	my $phylip = $self->_make_phylip( $taxa, @matrix );
	my $intree = $self->_make_intree( $taxa, $tree );
	my $binary = $self->_make_binary( $phylip );
	
	# execute
	my $string = $self->executable . $self->_setparams($binary,$intree);
	my $curdir = getcwd;
	chdir $self->work_dir;	
	system($string) and $self->warn("Couldn't run ExaML: $?");
	chdir $curdir;	
	return $self->_cleanup;
}

sub _cleanup {
	my $self = shift;
	my $dir  = $self->work_dir;
	my $out  = $self->outfile_name;
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
	return "${dir}/ExaML_result\.${out}";
}

sub _make_binary {
	my ( $self, $phylip ) = @_;
	my $binfile = File::Spec->catfile( $self->work_dir, $self->run_id . '-dat' );
	my ( $volume, $directories, $base ) = File::Spec->splitpath( $binfile );
	my $curdir = getcwd;
	chdir $self->work_dir;
	my @command = ( $self->parser, 
		'-m' => 'DNA', 
		'-s' => $phylip, 
		'-n' => $base,
		'>'  => File::Spec->devnull,		
		'2>' => File::Spec->devnull,
	);
	my $string = join ' ', @command;
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