#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';

# this is a unit test
BEGIN{ use_ok('Bio::Phylo::PhyLoTA::DAO'); }
my $obj = new_ok('Bio::Phylo::PhyLoTA::DAO');


