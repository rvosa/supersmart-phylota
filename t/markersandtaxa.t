#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';

# load the package
BEGIN { use_ok('Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa'); }

# create a new instance
my $mt = new_ok('Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa');

my $taxa = $Bin . '/testdata/species-allranks.tsv';

# parse a taxa table
my @records = $mt->parse_taxa_file($taxa);

ok (@records, 'successfully parsed taxa table');

# test subroutine get_root_taxon_level
my $rank = $mt->get_root_taxon_level(@records);
ok ($rank eq 'order', "root taxon is of rank $rank");

# test get_root_taxon_level with different species table
$taxa = $Bin . '/testdata/species.tsv';
@records = $mt->parse_taxa_file($taxa);
$rank = $mt->get_root_taxon_level(@records);
ok ($rank eq 'order', "root taxon is of rank $rank");


