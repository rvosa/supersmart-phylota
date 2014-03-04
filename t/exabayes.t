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

#ok( $exabayes->work_dir($wd), "work dir: $wd" );
#ok( $exabayes->executable( $config->EXABAYES_BIN ), "executable" );

ok( $exabayes->version == qv(1.2.1), "version" );
ok( $exabayes->m('DNA'), "model" );
ok( $exabayes->n('testrun'), "ruid" );
ok( $exabayes->s(1), "seed");
ok( $exabayes->w($wd), "working directory");
ok( $exabayes->C(2), "number of chains");
ok( $exabayes->R(2), "number of runs");

#$exabayes->run("$Bin/testdata.phy");
$exabayes->run("$Bin/testdata.xml");

print "Directory : $wd \n";
