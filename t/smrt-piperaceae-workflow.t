#!/usr/bin/perl

use warnings;

BEGIN {
	our @plan = $ENV{'TEST_SMRT'} ? 'no_plan' : 'skip_all' => 'env var TEST_SMRT not set';
}

use Test::More @plan;
use FindBin '$Bin';
use File::Temp qw(tempdir);
use File::Copy;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# This test runs the whole workflow for the piperaceae example and checks if a 
# final tree is produced. 

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new( '-level' => INFO, '-class' => 'main' );

# copy workflow and possible name- and fossil files into temporary directory
my $workdir = tempdir( CLEANUP => 1 );
my $srcdir = $ENV{'SUPERSMART_HOME'} . "/examples/piperaceae/";
chdir( $workdir );
my @files = glob("$srcdir/*");
copy( $_, $workdir ) or die "Copy failed: $!" for @files;

# run the whole workflow
$log->info("Running workflow from directory $srcdir. This may take a while.");
system ( "sh workflow.sh >/dev/null 2>&1" );

# check for final tree file
my $outfile = "$workdir/final.nex";

ok ( -e $outfile, "outfile exists" );
ok ( -s $outfile, "outfile not empty" );

# check if final tree can be read
my $tree = parse_tree(
	'-format' => 'nexus',
	'-file' => $outfile,
);

isa_ok ($tree, 'Bio::Phylo::Forest::Tree');
