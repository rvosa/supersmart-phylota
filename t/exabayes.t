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

ok( $exabayes->outfile_name('exabayes_out') );
ok( $exabayes->version == qv(1.2.1), "version" );
ok( $exabayes->m('DNA'), "model" );
ok( $exabayes->run_id('testrun'), "ruid" );
ok( $exabayes->s(1), "seed");
ok( $exabayes->work_dir( $wd ), "working directory" );
ok( $exabayes->C(1), "number of chains" );
ok( $exabayes->R(1), "number of runs" );
ok( $exabayes->parser("/usr/local/src/exabayes-1.2.1/bin/parser"), "set parser" );
ok( $exabayes->z(1), "quiet mode" );
ok( $exabayes->nodes(4), "number of nodes for parallel processing" );
ok( $exabayes->d(1), "performing dry runs");

# set MPI parameters
my $mpirun_bin = $config->MPIRUN_BIN;
ok( $exabayes->mpirun($mpirun_bin), "mpirun executable");
my $nodes      = $config->NODES;
ok( $exabayes->nodes($nodes), "number parallel nodes");

# set config file parameters
ok( $exabayes->numRuns(2), "number of Monte Carlo runs ");
ok( $exabayes->numCoupledChains(2), "coupled Monte Carlo chains ");
ok( $exabayes->numGen(2), "number of generations ");

# dry run very small example 
my $phylip = "$Bin/example.phy";
my $intree = "$Bin/example.dnd";
#ok( $exabayes->run('-phylip'=> $phylip, '-intree'=> $intree));

# dry run larger example
$phylip = "$Bin/../examples/primates/supermatrix.phy";
$intree = "$Bin/user.dnd";
$exabayes->run_id('testrun2');
ok( $exabayes->outfile_name('exabayes_out2'));
ok( $exabayes->run('-phylip'=> $phylip, '-intree'=> $intree));

# dry run from nexml file
$exabayes->n('testrun3');
ok( $exabayes->outfile_name('exabayes_out3'));
ok( $exabayes->run("$Bin/testdata.xml"));
