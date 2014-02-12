#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';

use Bio::Phylo::PhyLoTA::Service;

# this is a unit test
BEGIN{ use_ok('Bio::Phylo::PhyLoTA::DAO'); }

# test instanciating
my $schema = new_ok('Bio::Phylo::PhyLoTA::DAO');

# get ID for taxon name 
my $name = "Forsteronia refracta";
my $node = $schema -> resultset('Node')->find({ taxon_name => $name }); ##Hannes
my $id = $node->get_id();
is ($id, 387687, "test correct taxon id");
