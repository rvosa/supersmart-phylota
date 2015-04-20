package Bio::SUPERSMART::App::smrtutils::Command::OneZoomify;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse);

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: puts a 

=head1 NAME

OneZoomify - Makes a beadutiful OneZoom html figure from a phylogenetic tree

=head1 SYNOPSYS

smrt-utils onezoomify [-h ] [-v ] [-w <dir>] [-l <file>] [-y ] -t <file> [-o <file>] [-s ] 

=head1 DESCRIPTION

Creates a html file containing a OneZoom representation of a phylogenetic tree (see http://www.onezoom.org/).
The 'show' option can be used to open the figure immediately in the system's default internet browser.
This feature requires the perl package 'Browser::Open' to be installed.

=cut

sub options {
    
	my ($self, $opt, $args) = @_;		
	my $outfile_default = 'onezoom-tree.htm';
	return (
		['tree|t=s', 'file name of input tree (newick format)', { arg => 'file', mandatory => 1}],
		['outfile|o=s', "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => 'file'}],	
		['show|s', 'show generated figure in default browser, note that perl package Browser::Open must be installed for this option to work', {}],
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
	my $file = $ENV{SUPERSMART_HOME} . '/data/ONEZOOM/onezoom.htm';
	open my $fh, '<', $file;
	read $fh, my $oz, -s $fh;
	close $fh;
	
	my $placeholder = 'SUPERSMART_TREE_PLACEHOLDER';
	$oz =~ s/$placeholder/$treestr/g;

	open my $outfh, '>', $outfile or die $!;
	print $outfh $oz;
	close $outfh;
	
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
