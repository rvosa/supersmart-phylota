#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';
use Scalar::Util 'looks_like_number';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::TreeService;

# instantiate tree service
my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;

{
	# test nexus file
	my $file = "${Bin}/testdata/testconsense.nex";


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
}

{
	# test newick file
	my $file = "${Bin}/testdata/testconsense.dnd";


	# make consensus tree
	my $tree = $ts->consense_trees( '-infile' => $file, '-format' => 'newick' );
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
}