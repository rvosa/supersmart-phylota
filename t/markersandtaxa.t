#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';

# load the package
BEGIN { use_ok('Bio::SUPERSMART::Domain::MarkersAndTaxa'); }

# create a new instance
my $mt = new_ok('Bio::SUPERSMART::Domain::MarkersAndTaxa');

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

# test querying taxa table for taxon ids to get
#  ids for specified ranks

# get all species for subfamily "Cercopithecinae"
my $query = 9528;
my @ids = $mt->query_taxa_table( $query, "species", @records);
ok ( scalar(@ids) > 1, "found more than one species for subfamily Cercopithecinae in taxa table" );

my @query = (9605, 9596);
@ids = $mt->query_taxa_table( \@query, "species", @records);
ok ( scalar(@ids) > 1, "found more than one taxon for genera Homo and Pan in taxa table" );

# query for genera for a set of species
@query = (13515, 182256);
@ids = $mt->query_taxa_table( \@query, "genus", @records);

ok ( scalar(@ids) > 1, "found more than one genus for species 13515, 182256 in taxa table" );

@ids = $mt->query_taxa_table( \@query, "order", @records);
ok ( scalar(@ids) == 1, "found exactly one order for species 13515, 182256 in taxa table" );

# query all species and genera for order primates
my @ranks = ("species", "genus");
@ids = $mt->query_taxa_table( 9443, \@ranks, @records);
ok ( scalar(@ids) > 100, "found more than 100 species and genera for order primates" );

# test extracting subsets from fasta file
my $file = $Bin . '/testdata/aln.fa';

my %aln = $mt->parse_fasta_file( $file );
ok ( scalar(keys %aln) > 0, "could parse from fasta file" );

my @taxa = (446108, 290267, 446107);
my %subset = $mt->get_alignment_subset(\%aln, {'taxon' => \@taxa});

ok ( scalar(keys %subset) == 3, "could get subset of alignment");
