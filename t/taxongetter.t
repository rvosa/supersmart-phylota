#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Matrices::Matrix;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;

use Data::Dumper;

# BEWARE: this is a lengthy test to run!

my $config = Bio::Phylo::PhyLoTA::Config->new;

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new( '-level' => INFO, '-class' => 'main' );
$log->warn('BEWARE: this is a lengthy test to run!');

# the first tests: can we use and instantiate the MarkersAndTaxaSelector
BEGIN { use_ok('Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector'); }
my $mts = new_ok('Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector');
$log->VERBOSE(
	'-level' => WARN,
	'-class' => 'Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector',
);

SKIP: {
	skip "variable ROOT_TAXA_EXPANSION not set in phylota.ini" unless $config->ROOT_TAXA_EXPANSION;
	# get species for a given root taxon
	my @root_taxa_names = ("Gentianales", "Rubiaceae");

	my @taxa_names = $mts->expand_taxa(@root_taxa_names);
	cmp_ok (scalar @taxa_names, '>', 7000, "expand for root taxa");

	@root_taxa_names = ("Primates");
	@taxa_names = $mts->expand_taxa(@root_taxa_names);
	cmp_ok (scalar @taxa_names, '>', 300, "expand for root taxa");

	# 'Otolemur garnettii' has a subspecies, which should be 
	#  discarded if lowest level is "Species"
	@root_taxa_names = ("Otolemur garnettii");
	@taxa_names = $mts->expand_taxa(@root_taxa_names);
	# should be 'Otolemur garnettii' and NOT 'Otolemur garnettii garnettii'
	cmp_ok (scalar @taxa_names, '==', 1, "single taxon");
	cmp_ok ( $taxa_names[0], 'eq', "Otolemur garnettii", "test subspecies");
}

# we are going to read the text files in the results/specieslists dir,
my $file = $Bin . '/testdata/testnames.txt';
$log->info("going to read species names from $file");

# read names from file, clean line breaks
open my $fh, '<', $file or die $!;
my @names = <$fh>;
chomp(@names);
$log->info("read species names from $file");

# this will take some time to do the taxonomic name resolution in the
# database and with webservices
my @nodes = $mts->get_nodes_for_names(@names);
ok( @nodes, "found nodes" );

# test searching for taxon name that is not in database but 
# could be found using the webservice
my $taxon = "Alouatta fusca";
my @fusca_nodes = $mts->get_nodes_for_names($taxon);
# after TNRS searching webservice, Alouatta guariba should be found. 
is ($fusca_nodes[0]->taxon_name, "Alouatta guariba", "found species Alouatta guariba");

# search for taxon that cannot be found anywhere
my @no_nodes = $mts->get_nodes_for_names("xy");
is ( scalar @no_nodes, 0, "nonsense query yields no taxon");

# test getting the nodes for a whole taxon table 
my $mt =  Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
my $taxa = $Bin . '/testdata/species.tsv';
my @records = $mt->parse_taxa_file($taxa);
@nodes = $mts->get_nodes_for_table(@records);
ok( @nodes, "getting nodes for taxa table" );

# test getting nodes with a different taxon table
$taxa = $Bin . '/testdata/species-allranks.tsv';
@records = $mt->parse_taxa_file($taxa);
@nodes = $mts->get_nodes_for_table(@records);
ok( @nodes, "getting nodes for taxa table (many ranks)" );
