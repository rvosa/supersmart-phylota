package Bio::SUPERSMART::App::smrtutils::Command::Plot;
use strict;
use warnings;
use Template;
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

Creates a file with a graphical representation of a tree. Currently supported formats:

=over

=item Onezoom (html format, see L<http://www.onezoom.org>)

=item jsPhyloSVG (html format, see L<http://www.jsphylosvg.org>)

=item supersmart (html format)

=back

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
	my $width_default     = 800;
	my $height_default    = 600;
	my $style_default     = 'rectangular';
	my $bbmarkers_default = 'markers-backbone.tsv';
	my $taxa_default      = 'species.tsv';
	my $fossil_default    = 'fossils.tsv';
	my $clade_default     = 0;
	return (
		['tree|t=s', 'file name of input tree (newick format)', { arg => 'file', mandatory => 1 } ],
		['outfile|o=s', "name of the output file, defaults to '$outfile_default'", { default => $outfile_default, arg => 'file' } ],	
		['format|f=s', "output format, currently supported: OneZoom, jsPhyloSVG, SUPERSMART. Default: $format_default", { default => $format_default } ],
		['width|w=i', "visualization width. Default: $width_default", { default => $width_default } ],
		['height|h=i', "visualization height. Default: $height_default", { default => $height_default } ],
		['style|s=s', "visualization style (rectangular or circular). Default: $style_default", { default => $style_default } ],
		['markers|m=s', "backbone markers table. Default: $bbmarkers_default", { default => $bbmarkers_default } ],
		['taxa|i=s', "taxa table. Default: $taxa_default", { default => $taxa_default } ],
		['fossils|p=s', "fossil table. Default: $fossil_default", { default => $fossil_default } ],
		['clades|c=s', "Search cladeXXX folders in workdir. Default: $clade_default", { default => $clade_default } ],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		
	
	# we only have to check the 'infile' argument. 
	# if the infile is absent or empty, abort  
	my $file = $opt->tree;
	$self->usage_error('no tree argument given') if not $file;
	$self->usage_error('file $file does not exist') unless (-e $file);
	$self->usage_error('file $file is empty') unless (-s $file);			
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
# - show age ranges from *BEAST
# - show bootstrapped age ranges

sub run {
	my ($self, $opt, $args) = @_;
	my $outfile = $opt->outfile;
	my $logger  = $self->logger;	
	
	# read tree, clean labels
	my $tree = parse_tree(
		'-file'   => $opt->tree,
		'-format' => 'newick',
	);
	$tree->visit(sub{
		my $n = shift;
		my $name = $n->get_name;
		$name =~ s/_/ /g;
		$name =~ s/'(.+?)'/$1/;
		$n->set_name($name);
	});
	
	# compute coordinates and round to nearest integer
	my $drawer = Bio::Phylo::Treedrawer->new(
		'-width'  => $opt->width,
		'-height' => $opt->height,
		'-tree'   => $tree,
	);
	$drawer->compute_coordinates;
	$tree->visit(sub{
		my $n = shift;
		$n->set_x( int( $n->get_x + 0.5 ) );
		$n->set_y( int( $n->get_y + 0.5 ) );
	});
	
	# apply metadata for styling
	my $ds = Bio::Phylo::PhyLoTA::Service::DecorationService->new;
	$ds->apply_taxon_colors($opt->taxa,$tree) if -s $opt->taxa;
	$ds->apply_backbone_markers($opt->markers,$tree) if -s $opt->markers;
	$ds->apply_fossil_nodes($opt->fossils,$tree) if -s $opt->fossils;
	$ds->apply_clade_markers($opt->clades,$tree) if -d $opt->clades;
	
	# get template file
	my $tool = lc $opt->format;
	my $file = $ENV{'SUPERSMART_HOME'} . '/data/VISUALIZATION/' . $tool . '.tmpl';
		
	# instantiate template
	my $tt = Template->new( 'ABSOLUTE' => 1 );
	my $date = localtime();
	my %args = (
		'tree'    => $tree,
		'width'   => $opt->width,
		'height'  => $opt->height,
		'style'   => $opt->style,
		'date'    => $date,
		'command' => "$0 @ARGV",
	);
	my $oz = $tt->process( $file, \%args, $outfile ) or die $tt->error;
	$logger->info("DONE. Tree written to $outfile");
	return 1;
}

1;
