#!/usr/bin/perl
use strict;
use warnings;
use version;
use FindBin '$Bin';
use File::Temp 'tempdir';
use Test::More 'no_plan';
use Bio::Phylo::PhyLoTA::Config;

my $config = Bio::Phylo::PhyLoTA::Config->new;

my $wd = tempdir( 'CLEANUP' => 0 );

# test instantiation
BEGIN { use_ok('Bio::Tools::Run::Phylo::ExaBayes'); }
my $exabayes = new_ok('Bio::Tools::Run::Phylo::ExaBayes');

##ok( $exabayes->executable( $config->EXABAYES_BIN ), "executable" );

ok( $exabayes->version == qv(1.2.1), "version" );
ok( $exabayes->m('DNA'), "model" );
ok( $exabayes->n('testrun'), "ruid" );
ok( $exabayes->s(1), "seed");
ok( $exabayes->w($wd), "working directory");
ok( $exabayes->C(1), "number of chains");
ok( $exabayes->R(1), "number of runs");
ok( $exabayes->parser("/usr/local/src/exabayes-1.2.1/bin/parser"), "set parser");
ok( $exabayes->z(1), "quiet mode");
ok( $exabayes->d(1), "performing dry runs");

# dry run very small example 
my $phylip = "$Bin/example.phy";
my $intree = "$Bin/example.dnd";
ok( $exabayes->run('-phylip'=> $phylip, '-intree'=> $intree));

# dry run larger example
ok( $exabayes->run('-phylip'=> $phylip, '-intree'=> $intree));
$phylip = "$Bin/supermatrix.phy";
$intree = "$Bin/common.dnd";
ok( $exabayes->run('-phylip'=> $phylip, '-intree'=> $intree));

# dry run from nexml file
ok( $exabayes->run("$Bin/testdata.xml"));
