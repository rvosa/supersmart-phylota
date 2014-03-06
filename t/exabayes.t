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
ok( $exabayes->C(1), "number of chains");
ok( $exabayes->R(1), "number of runs");
ok( $exabayes->parser("/usr/local/src/exabayes-1.2.1/bin/parser"), "set parser");

my $phylip = "$Bin/supermatrix.phy";
my $intree = "$Bin/common.dnd";

#my $phylip = "$Bin/example.phy";
#my $intree = "$Bin/example.dnd";

my $parser = $exabayes->parser;
print "Parser = $parser \n";


$exabayes->run('-phylip'=> $phylip, '-intree'=> $intree);
##$exabayes->run("$Bin/testdata.xml");

print "Directory : $wd \n";
