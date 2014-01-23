#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

# process command line arguments
my $verbosity = WARN;
my @levels = qw[species genus family order class phylum kingdom];
my $root = 4055; # Gentianales
GetOptions(
	'root=i'   => \$root,
	'verbose+' => \$verbosity,
	'level=s'  => \@levels,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => [qw(
		main
		Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector
	)]
);
my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;

# print table header
print join("\t", 'name', @levels), "\n";

# fetch the root node object
my $root_node = $mts->find_node($root);
$log->info("fetched node object $root_node for ID $root");

# create hash of taxonomic levels so that when we walk up the
# taxonomy tree we can more easily check to see if we are at a
# level of interest
my %level = map { $_ => 1 } @levels;

# this will hold the higher level taxa for all nodes we visit
my %levels_for;

# build the path up to the root
my $node      = $root_node;
my $root_rank = $root_node->rank;
my @rootpath;
while( $node ) {
	my $rank = $node->rank;
	my $id   = $node->ti;
	if ( $level{$rank} && $id != $root ) {
		
		# because we are now traversing from child to parent we need to push
		# into the path array (i.e. append to it), whereas later on when we
		# traverse from parent to child we need to unshift (i.e. pre-pend)
		push @rootpath, $id; # NB!
	}
	$node = $node->get_parent;
}
$log->info("fetched path from $root to root: @rootpath");

# do a depth-first, pre-order traversal, growing the node paths
# recursively, and writing them out when we arrive at a tip
$root_node->visit_depth_first(
	'-pre' => sub {
		my $node = shift;
		my $id   = $node->ti;
		my $rank = $node->rank;
		$log->info("focal node $id has rank $rank");
		
		# this happens once the recursion starts to attempt to
		# process sister clades of the root clade, which it shouldn't
		if ( $rank eq $root_rank && $id != $root ) {
			$log->info("arrived at sister-of-root-clode $id, quitting...");
			exit(0);
		}
		
		# fetch the path from the parent, if any
		my @path;
		if ( ( my $parent = $node->get_parent ) && $id != $root ) {
			@path = @{ $levels_for{$parent->ti} };
			$log->info("fetched path @path from parent");
		}
		
		# extend the path if we are at a level of interest
		if ( $level{$rank} ) {
			
			# here we pre-pend, see above
			unshift @path, $id; # NB!
			$log->info("going to store rank '$rank' as a node in the path");
		}
		
		# print the path if
		if ( $rank eq $levels[0] ) {
			print join("\t", $node->taxon_name, @path, @rootpath), "\n";
		}
		$levels_for{$id} = \@path;
	},
	'-post' => sub {
		my $node = shift;
		my $id   = $node->ti;
		if ( $id == $root ) {
			$log->info("returned to root $root in post-order, quitting...");
			exit(0);
		}
	},
);