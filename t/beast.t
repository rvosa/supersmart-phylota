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
ok( $beast->collapse_species(1), "collapse species");
ok( $beast->rebuild(1), "rebuild input files");

my ($handle, $file) = tempfile( 'OPEN' => 1 );
close $handle;

ok( $beast->outfile_name( $file ), "outfile");
ok( $beast->run( "$Bin/testdata/testdata.xml" ), "run BEAST");


# Caution, lengthy test ! 
# second example, test if subspecies, varietas and forma 
# will be present in output tree
my $runlengthy = 0;
SKIP : {
	skip "running second beast example would take too long " unless $runlengthy;
	ok ($beast->collapse_species(0) == 0, "do not collapse species");
	ok( $beast->chain_length(10) == 10, "chain length" );
	ok( $beast->sample_freq(10) == 10, "sample_freq" );
	ok( $beast->overwrite(1), "overwrite");
	ok( $beast->run( "$Bin/testdata/testdata_cucumis.xml" ), "run BEAST");
};
