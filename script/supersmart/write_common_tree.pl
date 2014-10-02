#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO qw(parse_tree unparse);
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;

=head1 NAME

write_common_tree.pl - writes polytomous taxonomy tree

=head1 SYNOPSYS

 $ write_common_tree.pl --infile=<taxa table> > <outfile>

=head1 DESCRIPTION

Given a table of reconciled taxa, writes the 'common tree' that connects these taxa
in the underlying taxonomy. The resulting tree description is written to STDOUT. By
default this is in Newick syntax, with labels on interior nodes (including all available
taxonomic levels up to the root), and with the taxon IDs as terminal labels.

=cut

# process command line arguments
my $infile     = '-'; # read from STDIN by default
my $outformat  = 'newick'; # write newick by default
my $nodelabels = 1; # write internal node labels by default
my $verbosity  = INFO; # low verbosity
GetOptions(
	'infile=s'    => \$infile,
	'outformat=s' => \$outformat,
	'verbose+'    => \$verbosity,
	'nodelabels'  => $nodelabels,
);

# instantiate helper objects
my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
my $mt =  Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;

my $log = Bio::Phylo::Util::Logger->new(
	'-class' => [ 'main', 'Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector' ],
	'-level' => $verbosity
);

# parse the taxa file 
my @taxatable = $mt->parse_taxa_file($infile);

# instantiate nodes from infile
my @nodes = $mts->get_nodes_for_table( @taxatable );

# compute common tree
my $tree = $mts->get_tree_for_nodes(@nodes);
$log->debug("done computing common tree");

my $newick = $tree->to_newick;

# create node labels
$tree->visit(sub{
	my $node = shift;
	if ( $nodelabels or $node->is_terminal ) {
		my $label = $node->get_guid;		
		$node->set_name( $label );
	}
});

# write output
print unparse(
	'-format'     => $outformat,
	'-phylo'      => $tree,
	'-nodelabels' => $nodelabels,
);