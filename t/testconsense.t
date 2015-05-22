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
	# test newick file
	my $file = "${Bin}/testdata/testconsense.dnd";


	# make consensus tree
	my $tree = $ts->consense_trees( '-infile' => $file, '-format' => 'newick' );
	ok( $tree, "have tree $tree" );

	# check if the semantic annotations are numerically correct
	my $tip = $tree->get_by_name('taxon_4');
	
	my %exp = (		
		'height_95_HPD_min'  => 0.0,
		'height_95_HPD_max'  => 1.7763568394002505E-15,
		'length_range_min'   => 0.5128205002832722,
		'length_range_max'   => 0.9142517343483778,
		'height_median'      => 0.0,
		'length_95_HPD_min'  => 0.5128205002832722,
		'length_95_HPD_max'  => 0.9142517343483778,
		'height'             => 1.9737298215558337E-16,
		'height_range_min'   => 0.0,
		'height_range_max'   => 1.7763568394002505E-15,
		'length'             => 0.6526713962566253,
		'length_median'      => 0.6258248261274582,
	);
	for my $key ( keys %exp ) {
		my $obs = $tip->get_meta_object( 'fig:' . $key );
		ok( $obs == $exp{$key}, "$key $obs == $exp{$key}" );
	}
}