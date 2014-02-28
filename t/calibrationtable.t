#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use FindBin '$Bin';
use Bio::Phylo::IO 'parse_tree';

my $ctc_class = 'Bio::Phylo::PhyLoTA::Service::CalibrationTableCreator';
my $fdg_class = 'Bio::Phylo::PhyLoTA::Service::FossilDataGetter';
my $ct_class = 'Bio::Phylo::PhyLoTA::Domain::CalibrationTable';

# test to see if the service and data classes can be loaded dynamically
use_ok($ctc_class);
use_ok($fdg_class);
use_ok($ct_class);

# test to see if we can instantiate objects of the service class
my $ctc = new_ok($ctc_class);
my $fdg = new_ok($fdg_class);

# fetch the associated config object
my $config = $ctc->config;
isa_ok( $config, 'Bio::Phylo::PhyLoTA::Config' );

# test to see if a fossil table file can be found
my $file = "$Bin/../examples/gentianales/fossils3.tsv";
ok( -e $file, "$file exists" );

# read the fossil table
my @fossils = $fdg->read_fossil_table( $file );
ok( scalar(@fossils), "read ".scalar(@fossils)." fossils" );

# find calibration points
my @points;
for my $f ( @fossils ) {    
    if ( my $n = $ctc->find_calibration_point($f) ) {
        isa_ok( $n, 'Bio::Phylo::PhyLoTA::Domain::FossilData' );	
        push @points, $n;
    }
}
ok( scalar(@points), "converted to ".scalar(@points)." points" );

# get a tree to test creating a calibration table (Gentianales tree)
my $tree = parse_tree( 
    '-format'     => 'newick', 
    '-file'       => "$Bin/testtree.dnd", 
    '-as_project' => 1 ,
    '-keep-whitespace' => 1
);
isa_ok($tree, 'Bio::Phylo::Forest::Tree');

# make calibration table from fossils
$config->fossil_best_practice_cutoff(0.0);
my $table = $ctc->create_calibration_table( $tree, @points );
isa_ok($table,  'Bio::Phylo::PhyLoTA::Domain::CalibrationTable');
my @rows = $table->get_rows();

# set higher best practice score cutoff and see if there is
# less rows generated in a calibration table
$config->fossil_best_practice_cutoff(3.0);
my $table2 = $ctc->create_calibration_table( $tree, @points );
my @rows2 = $table2->get_rows();

ok(scalar(@rows) > scalar(@rows2), 'Less rows in calibration table when cutoff is increased');

