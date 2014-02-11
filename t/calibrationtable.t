#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use FindBin '$Bin';

my $ctc_class = 'Bio::Phylo::PhyLoTA::Service::CalibrationTableCreator';
my $fdg_class = 'Bio::Phylo::PhyLoTA::Service::FossilDataGetter';

# test to see if the service classes can be loaded dynamically
use_ok($ctc_class);
use_ok($fdg_class);

# test to see if we can instantiate objects of the service class
my $ctc = new_ok($ctc_class);
my $fdg = new_ok($fdg_class);

# fetch the associated config object
my $config = $ctc->config;
isa_ok( $config, 'Bio::Phylo::PhyLoTA::Config' );

# test to see if a fossil table file can be found
my $file = "$Bin/../examples/gentianales/fossils2.tsv";
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