package Bio::SUPERSMART::App::smrtutils::Command::Plot;
use strict;
use warnings;
use Template;
use Bio::Phylo::IO qw(parse_tree);
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

=item Onezoom (html format, see www.onezoom.org)

=back

=cut

my %formats = ( 
	'onezoom' => 'html', 
);

sub options {
    
	my ($self, $opt, $args) = @_;
	my $format_default  = 'onezoom';
	my $outfile_default = "$format_default.$formats{$format_default}";
	my $width_default   = 800;
	my $height_default  = 600;
	my $style_default   = 'rectangular';
	return (
		['tree|t=s', 'file name of input tree (newick format)', { arg => 'file', mandatory => 1 } ],
		['outfile|o=s', "name of the output file, defaults to '$outfile_default'", { default => $outfile_default, arg => 'file' } ],	
		['format|f=s', "output format, currently supported: OneZoom, jsPhyloSVG. Default: $format_default", { default => $format_default } ],
		['width|w=i',"visualization width. Default: $width_default", { default => $width_default } ],
		['height|h=i',"visualization height. Default: $height_default", { default => $height_default } ],
		['style|s=s',"visualization style (rectangular or circular). Default: $style_default", { default => $style_default } ],
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

sub run {
	my ($self, $opt, $args) = @_;
	my $outfile = $opt->outfile;
	my $logger  = $self->logger;
	
	# read tree
	my $tree = parse_tree(
		'-file'   => $opt->tree,
		'-format' => 'newick',
	);
	
	# get template file and add newick string
	my $tool = lc $opt->format;
	my $file = $ENV{'SUPERSMART_HOME'} . '/data/VISUALIZATION/' . $tool . '.tmpl';
		
	# instantiate template
	my $tt = Template->new( 'ABSOLUTE' => 1 );	
	my %args = (
		'tree'   => $tree,
		'width'  => $opt->width,
		'height' => $opt->height,
		'style'  => $opt->style,
	);
	my $oz = $tt->process( $file, \%args, $outfile ) or die $tt->error;
	$logger->info("DONE. Tree written to $outfile");
	return 1;
}

1;
