#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';
use Scalar::Util 'looks_like_number';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::TreeService;

# configure helpers and files
my $file = "${Bin}/testdata/testconsense.nex";
my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;

# make consensus tree
my $tree = $ts->consense_trees( '-infile' => $file );
ok( $tree, "have tree $tree" );

# verify annotations
$tree->visit(sub{
	my $n = shift;
	my $ann = $n->get_generic('treeannotator');
	ok( $ann, "$ann exists" );
	ok( ref $ann eq 'HASH', "$ann is a hash" );
	for my $key ( keys %$ann ) {
		my $val = $ann->{$key};
		if ( ref $val ) {
			ok( ref $val eq 'ARRAY', "$val is an array" );
			for my $v ( @$val ) {
				ok( looks_like_number $v, "$v is a number" );
			}
		}
		else {
			ok( looks_like_number $val, "$val is a number" )
		}
	}

});