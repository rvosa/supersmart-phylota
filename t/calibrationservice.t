#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';

BEGIN {
	use_ok('Bio::Phylo::PhyLoTA::Service::CalibrationService');
}

my $cs = new_ok('Bio::Phylo::PhyLoTA::Service::CalibrationService');

# read a PROTEUS compatible tab-separated spreadsheet with fossils,
# returns an array of Bio::Phylo::PhyLoTA::Domain::FossilData objects
my @fossils = $cs->read_fossil_table("$Bin/../examples/primates/fossils.tsv");
for my $f ( @fossils ) {
	isa_ok( $f, 'Bio::Phylo::PhyLoTA::Domain::FossilData' );
}

# identify the calibration points the nodes should attach to
my @identified = map { $cs->find_calibration_point($_) } @fossils;
for my $i ( @identified ) {
	isa_ok( $i, 'Bio::Phylo::PhyLoTA::Domain::FossilData' );
}

# read a newick tree
my $tree = parse(
	'-format' => 'newick',
	'-file'   => 'tree.dnd',
	'-as_project' => 1,
);
isa_ok( $tree, 'Bio::Tree::TreeI' );

# creates a Bio::Phylo::PhyLoTA::Domain::CalibrationTable
my $ct = $cs->create_calibration_table( $tree, @identified );
isa_ok( $ct, 'Bio::Phylo::PhyLoTA::Domain::CalibrationTable' );

# calibrate the tree
# my $chronogram = $cs->calibrate_tree($tree,$ct);
