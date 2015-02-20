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

# search for smilodon
$name = 'Smilodon populator';
@nodes = $service->search_node( { taxon_name => $name } )->all;
ok ( scalar @nodes > 0, "found at least one node for $name");

# search for Echinochloa crus-galli
$name = "Echinochloa crus-galli";
@nodes = $service->search_node( { taxon_name => $name } )->all;
ok ( scalar @nodes > 0, "found at least one node for $name");

# search for Bletia sp. Chase O-581
$name = "Bletia sp. Chase O-581";
@nodes = $service->search_node( { taxon_name => $name } )->all;
ok ( scalar @nodes > 0, "found at least one node for $name");

# search for x Brassolaeliocattleya 'Sung Ya Green'
$name = "x Brassolaeliocattleya 'Sung Ya Green'";
@nodes = $service->search_node( { taxon_name => $name } )->all;
ok ( scalar @nodes > 0, "found at least one node for $name");

# search for Chagunius chagunio
#$name = "Chagunius chagunio";
#@nodes = $service->search_node( { taxon_name => $name } )->all;
#ok ( scalar @nodes > 0, "found at least one node for $name");

# search for a genus: Gobius
$name = "Gobius";
@nodes = $service->search_node( { taxon_name => $name } )->all;
ok ( scalar @nodes > 0, "found at least one node for $name");


