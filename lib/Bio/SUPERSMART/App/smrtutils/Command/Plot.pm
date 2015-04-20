package Bio::SUPERSMART::App::smrtutils::Command::Plot;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse);

use Template;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: puts a 

=head1 NAME

Plot - Visualizes a tree

=head1 SYNOPSYS



=head1 DESCRIPTION

Creates a file with a graphical representation of a tree. Currently supported formats:

=over

=item Onezoom (html format, see www.onezoom.org)

=back

=cut

my %formats = ('onezoom' => 'htm');

sub options {
    
	my ($self, $opt, $args) = @_;
	my $format_default = 'onezoom';
	my $outfile_default = "$format_default.$formats{$format_default}";
	return (
		['tree|t=s', 'file name of input tree (newick format)', { arg => 'file', mandatory => 1}],
		['outfile|o=s', "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => 'file'}],	
		['format|f=s', "output format, currently supported: OneZoom. Default: $format_default", {default=>$format_default} ],
		['show|s', 'show generated figure in default tool, note that perl package Browser::Open must be installed for this option to work with OneZoom', {}],
	    );	
}

sub validate {
	my ($self, $opt, $args) = @_;		
	
	# We only have to check the 'infile' argument. 
	#  If the infile is absent or empty, abort  
	my $file = $opt->tree;
	$self->usage_error('no tree argument given') if not $file;
	$self->usage_error('file $file does not exist') unless (-e $file);
	$self->usage_error('file $file is empty') unless (-s $file);			
}

sub run {
	my ($self, $opt, $args) = @_;    
	my $treefile = $opt->tree;
	my $outfile = $opt->outfile;
	my $logger = $self->logger;
	
	# read tree
	my $treestr = parse(
		'-file'   => $treefile,
		'-format' => 'newick',
	    )->first->to_newick;
	
	# get onezoom template file and add newick string
	my $tool = $opt->format;
	my $file = $ENV{SUPERSMART_HOME} . '/data/VISUALIZATION/' . $tool . '.tmpl';
	
	
	my %args = (tree=>"$treestr");
	my $tt = Template->new(ABSOLUTE=>1);
	my $oz = $tt->process($file, \%args, $outfile) or die $tt->error;
	
	if ( $opt->show ) {
		
		# make Browser::Open a 'soft' dependency; it most likely won't work on VMs and cluster computers 
		my $bo = eval {
			require Browser::Open;
			1;
		};
		if ( $bo ) {
			Browser::Open::open_browser( $outfile ); 
		}
		else {
			$logger->warn('Browser::Open must be installed in order to show OneZoom tree in browser.')
		}
	}

	$logger->info("DONE. Onezoom tree written to $outfile");

	return 1;
}

1;
