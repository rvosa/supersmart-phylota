#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';

BEGIN{ use_ok('Bio::Phylo::PhyLoTA::Service'); }

my $service = new_ok ('Bio::Phylo::PhyLoTA::Service');

# test retreiving Trichosanthes cucumerina from database
my $name = 'Trichosanthes cucumerina';
my @nodes = $service->search_node( { taxon_name => $name } )->all;
 
my $node = shift @nodes;
my $id = $node->get_id();
is ( $id, 50543, "found $name");

# search for some taxon names
my @names = ("Smilodon populator", "Echinochloa crus-galli", 
			"Bletia sp. Chase O-581", "x Brassolaeliocattleya 'Sung Ya Green'",
			"Chagunius chagunio", "Schismatorhynchos nukta", "Garra waterloti",
			"Gobius");

for my $name (@names) {
	@nodes = $service->search_node( { taxon_name => $name } )->all;
	ok ( scalar @nodes > 0, "found at least one node for $name");
}
