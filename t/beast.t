#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';
use File::Temp qw(tempfile);
use Bio::Phylo::PhyLoTA::Config;

my $config = Bio::Phylo::PhyLoTA::Config->new;

# test instantiation
BEGIN { use_ok('Bio::Tools::Run::Phylo::StarBEAST'); }
my $beast = new_ok('Bio::Tools::Run::Phylo::StarBEAST');

ok( $beast->executable( $config->BEAST_BIN ), "executable" );
ok( $beast->chain_length(100) == 100, "chain length" );
ok( $beast->sample_freq(100) == 100, "sample_freq" );
ok( $beast->overwrite(1), "overwrite");

my (undef, $file) = tempfile(OPEN=>0);

ok( $beast->outfile_name( $file ), "outfile");
ok( $beast->run( "$Bin/testdata/testdata.xml" ), "run BEAST");