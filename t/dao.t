#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use Bio::SUPERSMART::Service;

# this is a unit test
BEGIN{ use_ok('Bio::Phylo::PhyLoTA::DAO'); }

# test instantiating
my $schema = new_ok('Bio::Phylo::PhyLoTA::DAO');

# get ID for taxon name 
my $name = "Forsteronia refracta";
my $node = $schema->resultset('Node')->find({ taxon_name => $name }); 
isa_ok($node, 'Bio::Phylo::PhyLoTA::DAO::Result::Node');
my $id = $node->get_id();
is($id, 387687, "test correct taxon id");

# Taxon Periploca has more than one row in the nodes table
$name = "Periploca";
my @nodes = $schema->resultset('Node')->search({ taxon_name => $name })->all; 
is (scalar @nodes, 2, "test retrieving multiple rows");

for my $n ( @nodes ){
    isa_ok( $n, 'Bio::Phylo::PhyLoTA::DAO::Result::Node');
}

# try retrieving Trichosanthes cucumerina. This should be in the database;
$name = "Trichosanthes cucumerina";
@nodes = $schema->resultset('Node')->search({ taxon_name => $name }); 
$node = shift @nodes;
$id = $node->get_id();
is ( $id, 50543, "found $name");
