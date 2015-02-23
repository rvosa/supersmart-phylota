#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

my @names = (
	'Homo sapiens',
	'Pan paniscus',
	'Pan troglodytes',
	'Gorilla gorilla',
	'Pongo pygmaeus',
);

my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
my @nodes = $mts->get_nodes_for_names(@names);
my $tree  = $mts->get_tree_for_nodes(@nodes);
ok( my $newick = $tree->to_newick( '-nodelabels' => 1 ) );
