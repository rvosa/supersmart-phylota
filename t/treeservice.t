#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa 'parse_taxa_file';

# load the package
BEGIN { use_ok('Bio::Phylo::PhyLoTA::Service::TreeService'); }

# create a new instance
my $ts = new_ok('Bio::Phylo::PhyLoTA::Service::TreeService');

my $tf = $Bin . '/testdata/testtree.dnd';
my $tree = parse_tree(
	'-file'   => $tf,
	'-format' => 'newick',
	'-as_project' => 1,
    );

my $taxa = $Bin . '/testdata/testspecieslist.tsv';

# parse taxon mapping
my $mt = 'Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa';
my @records = $mt->parse_taxa_file($taxa);

# test counting monophyletic families
#$ts->_count_monophyletic_groups($tree, @records);

# test rerooting a tree
#my $rerooted = $ts->reroot_tree($tree, @records);

#isa_ok ($rerooted, 'Bio::Tree::TreeI');

# rerooted tree should be different
#ok(! ($tree->to_newick eq $rerooted->to_newick));

# parse phylip files for tree tips
my $phylip_small = "$Bin/testdata/testmatrix-small.phy";
my @names_small = $ts->read_tipnames($phylip_small);
is(scalar @names_small, 5, "read_tipnames from phylip file");

my $phylip_large = "$Bin/testdata/testmatrix.phy";
my @names_large = $ts->read_tipnames($phylip_large);
is($names_large[0], "Orthogeomys_heterodus", "read_tipnames from phylip file");

# build a consensus tree
my $treefile  = "$Bin/testdata/testclade.nex";
my $consensus = $ts->consense_trees( '-infile' => $treefile );
isa_ok($consensus,'Bio::Tree::TreeI');

# parse the backbone tree
my $bbfile = "$Bin/testdata/testtree-labelled.dnd";
my $bbtree = parse_tree(
	'-file'   => $bbfile,
	'-format' => 'newick',
	'-as_project' => 1,
);
my $grafted = $ts->graft_tree( $bbtree, $consensus );
isa_ok($grafted,'Bio::Tree::TreeI');
