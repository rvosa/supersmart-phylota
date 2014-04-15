#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa 'parse_taxa_file';

my $primates = "$Bin/../examples/primates";

# load the package
BEGIN { use_ok('Bio::Phylo::PhyLoTA::Service::TreeService'); }

# create a new instance
my $ts = new_ok('Bio::Phylo::PhyLoTA::Service::TreeService');

my $tf = "$primates/megatree.dnd";
my $tree = parse_tree(
	'-file'   => $tf,
	'-format' => 'newick',
	'-as_project' => 1,
    );

my $taxa = "$primates/species.tsv";

# parse taxon mapping
my $mt = 'Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa';
my @records = $mt->parse_taxa_file($taxa);

my $rerooted = $ts->reroot_tree($tree, @records);

isa_ok ($rerooted, 'Bio::Tree::TreeI');

# rerooted tree should be different
ok(! ($tree->to_newick eq $rerooted->to_newick));

# parse phylip files for tree tips
my $phylip_small = "$Bin/example.phy";
my @names_small = $ts->read_tipnames($phylip_small);
is(scalar @names_small, 5, "read_tipnames from phylip file");

my $phylip_large = "$Bin/example2.phy";
my @names_large = $ts->read_tipnames($phylip_large);
is($names_large[0], "Orthogeomys_heterodus", "read_tipnames from phylip file");

# build a consensus tree
my $treefile  = "$primates/clade3/clade3.nex";
my $consensus = $ts->consense_trees( '-infile' => $treefile );
isa_ok($consensus,'Bio::Tree::TreeI');

# parse the backbone tree
my $bbfile = "$primates/supermatrix-cover1-rerooted-remapped-calibrated.dnd";
my $bbtree = parse_tree(
	'-file'   => $bbfile,
	'-format' => 'newick',
	'-as_project' => 1,
);
my $grafted = $ts->graft_tree( $bbtree, $consensus );
isa_ok($grafted,'Bio::Tree::TreeI');
