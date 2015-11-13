package Bio::SUPERSMART::App::smrtutils::Command::Plot;
use strict;
use warnings;
use Template;
use File::Spec;
use List::Util 'sum';
use List::MoreUtils qw(pairwise uniq all);
use Bio::Phylo::Treedrawer;
use Bio::Phylo::IO qw(parse_tree);
use Bio::Phylo::PhyLoTA::Service::DecorationService;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);


# ABSTRACT: plots a phylogenetic tree

=head1 NAME

Plot - Visualizes a tree

=head1 SYNOPSYS

=head1 DESCRIPTION

Creates a file with a graphical representation of a tree. 

Currently supported formats:

- Onezoom (html format, see L<http://www.onezoom.org>)
- jsPhyloSVG (html format, see L<http://www.jsphylosvg.org>)
- supersmart (html format)

=cut

my %formats = ( 
	'onezoom'    => 'html', 
	'supersmart' => 'html',
	'jsphylosvg' => 'html',
);

sub options {
	my ($self, $opt, $args) = @_;
	my $format_default    = 'supersmart';
	my $outfile_default   = "$format_default.$formats{$format_default}";
	my $width_default     = 1200;
	my $style_default     = 'rectangular';
	my $bbmarkers_default = 'markers-backbone.tsv';
	my $taxa_default      = 'species.tsv';
	my $fossil_default    = 'fossils.tsv';
	my $tree_default      = 'final.nex';
	my $clade_default     = 'markers-clades.tsv';
	return (
		['tree|t=s', "file name of input tree (NEXUS/figtree format, default is '$tree_default')", { default => $tree_default, arg => 'file' } ],
		['outfile|o=s', "name of the output file, defaults to '$outfile_default'", { default => $outfile_default, arg => 'file' } ],	
		['format|f=s', "output format, currently supported: OneZoom, jsPhyloSVG, SUPERSMART. Default: '$format_default'", { default => $format_default, arg => 'string' } ],
		['width|w=i', "visualization width. Default: $width_default", { default => $width_default, arg => 'int' } ],
		['height|e=i', "visualization height", { arg => 'int' } ],
		['style|s=s', "visualization style (rectangular or circular). Default: $style_default", { default => $style_default, arg => 'file' } ],
		['markers|m=s', "backbone markers table. Default: $bbmarkers_default", { default => $bbmarkers_default, arg => 'file' } ],
		['taxa|i=s', "taxa table. Default: $taxa_default", { default => $taxa_default, arg => 'file' } ],
		['fossils|p=s', "fossil table. Default: $fossil_default", { default => $fossil_default, arg => 'file' } ],
		['clades|c=s', "clade markers table. Default: $clade_default", { default => $clade_default, arg => 'file' } ],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		
	
	# we have to check the 'infile' argument. 
	# if the infile is absent or empty, abort  
	my $file = $opt->tree;
	$self->usage_error("tree file $file does not exist") unless -e $file;
	$self->usage_error("tree file $file is empty") unless -s $file;			

	# needs to be positive
	if ( $opt->width < 0 ) {
		$self->usage_error("width must be a positive integer, not ".$opt->width);
	}
	if ( defined( $opt->height ) and $opt->height < 0 ) {
		$self->usage_error("height must be a positive integer, not ".$opt->height);
	}

	# check if exists
	if ( not -e $opt->markers ) {
		$self->logger->info("backbone marker table not found, will not process ".$opt->markers); 
	}
	if ( not -e $opt->taxa ) {
		$self->logger->info("taxa table not found, will not process ".$opt->taxa);
	}
	if ( not -e $opt->fossils ) {
		$self->logger->info("fossil table not found, will not process ".$opt->fossils);
	}
	if ( not -e $opt->clades ) {
		$self->logger->info("clade marker table not found, will not process ".$opt->clades);
	}
}

# here is my wish list for this:
# + show which taxa were exemplars (bold taxon name)
# + show the backbone branches (dashed)
# + show bootstraps as fraction, e.g. 85/100
# + show posteriors as probability 
# + show higher taxa as evolving branch colors
# + mark fossil nodes (maybe mark ranges?)
# + show branch coverage as line thickness
# + make taxa clickable (popup with link to NCBI taxonomy and to used sequences)
# + show age ranges from *BEAST
# + show bootstrapped age ranges

sub run {
	my ($self, $opt, $args) = @_;
	my $outfile = $opt->outfile;
	my $logger  = $self->logger;	
	my $ntax    = 0;
	
	# read tree, clean labels
	$logger->info("Going to read tree file ".$opt->tree);
	my $tree = parse_tree(
		'-file'   => $opt->tree,
		'-format' => 'figtree',
	);
	$tree->visit(sub{
		my $n = shift;
		my $name = $n->get_name;
		$name =~ s/_/ /g;
		$name =~ s/'(.+?)'/$1/;
		$n->set_name($name);
		$ntax++ if $n->is_terminal;
	});
	
	# compute coordinates and round to nearest integer
	my $height = $opt->height || $ntax * 25;
	$logger->info("Going to compute node coordinates");
	my $drawer = Bio::Phylo::Treedrawer->new(
		'-width'  => $opt->width,
		'-height' => $height,
		'-tree'   => $tree,
	);
	$drawer->compute_coordinates;
	my ( @present, @pptu );
	$tree->visit(sub{
		my $n = shift;
		my $x = $n->get_x;
		push @present, $x if $n->is_terminal;
		my $l = $n->get_branch_length || 0;
		if ( $l > 0 and $n->get_parent ) {
			push @pptu, ( $x - $n->get_parent->get_x ) / $l;
		}
		$n->set_x( int( $x + 0.5 ) );
		$n->set_y( int( $n->get_y + 0.5 ) );
	});
	
	# apply metadata for styling
	my $ds = Bio::Phylo::PhyLoTA::Service::DecorationService->new;
	if ( -e $opt->taxa and -s $opt->taxa ) {
		$logger->info("Going to apply taxon colors from file ".$opt->taxa);
		$ds->apply_taxon_colors($opt->taxa,$tree);
	}
	if ( -e $opt->markers and -s $opt->markers ) {
		$logger->info("Going to apply backbone markers from file ".$opt->markers);
		$ds->apply_markers($opt->markers,$tree,'backbone');
	}
        if ( -e $opt->clades and -s $opt->clades ) {
                $logger->info("Going to apply clade markers from file ".$opt->clades);
                $ds->apply_markers($opt->clades,$tree,'clade');
        }
	if ( -e $opt->fossils and -s $opt->fossils ) {
		$logger->info("Going to apply fossil annotations from file ".$opt->fossils);
		$ds->apply_fossil_nodes($opt->fossils,$tree);
	}
	
	# get template file
	my $tool = lc $opt->format;
	my $file = $ENV{'SUPERSMART_HOME'} . '/data/VISUALIZATION/' . $tool . '.tmpl';
		
	# instantiate template
	my $tt = Template->new( 
		'ABSOLUTE'     => 1, 
		'INCLUDE_PATH' => $ENV{'SUPERSMART_HOME'} . '/data/VISUALIZATION/' 
	);
	my $date = localtime();
	my ( $vol, $dirs, $script ) = File::Spec->splitpath( $0 );	
	my %args = (
		'tree'    => $tree,
		'width'   => $opt->width,
		'height'  => $height,
		'style'   => $opt->style,
		'date'    => $date,
		'command' => "$script @ARGV",
		'present_x_coord'   => ( sum(@present)/scalar(@present) ),
		'pix_per_time_unit' => ( sum(@pptu)/scalar(@pptu) ),
	);
	my $oz = $tt->process( $file, \%args, $outfile ) or die $tt->error;
	$logger->info("DONE. Tree written to $outfile");
	return 1;
}

1;
