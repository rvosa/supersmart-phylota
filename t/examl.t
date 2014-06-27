#!/usr/bin/perl
use strict;
use warnings;
use version;
use FindBin '$Bin';
use File::Temp 'tempdir';
use Test::More 'no_plan';
use Bio::Phylo::PhyLoTA::Config;

my $wd = tempdir( 'CLEANUP' => 1 );
my $config = Bio::Phylo::PhyLoTA::Config->new;

# test instantiation
BEGIN { use_ok('Bio::Tools::Run::Phylo::ExaML'); }
my $examl = new_ok('Bio::Tools::Run::Phylo::ExaML');

ok( $examl->outfile_name('out'), "out file" );
ok( $examl->work_dir($wd), "work dir: $wd" );
ok( $examl->executable( $config->EXAML_BIN ), "executable" );
ok( $examl->version == qv(1.0.0), "version" );
ok( $examl->parser( $config->EXAML_PARSER_BIN ), "parser" );
ok( $examl->m('GAMMA'), "model" );
ok( $examl->quiet(1), "set quiet" );
ok( -e $examl->run( "$Bin/testdata.xml" ) );

