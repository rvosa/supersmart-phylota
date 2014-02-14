#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use Bio::Phylo::Matrices::Matrix;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;

# instantiate helper objects
my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
my $conf = Bio::Phylo::PhyLoTA::Config->new;

# this fetches the smallest containing cluster around the seed
# sequence with gi 326632174, which is a lemur cytb sequence
my @sequences = $sg->get_smallest_cluster_for_sequence(326632174);
my @filtered  = $sg->filter_seq_set(@sequences);

# if you are brave, you can use all of these, which have to be in the path:
# mafft clustalw kalign muscle probalign probcons tcoffee amap
for my $tool ( qw(muscle) ) {
	$conf->MSA_TOOL($tool); # runtime configuration change

	# align using whatever aligner is specified
	my $align = $sg->align_sequences(@filtered);

	# converts a Bio::AlignI object (as returned by muscle) into
	# a matrix object that can be written to nexus
	my $matrix = Bio::Phylo::Matrices::Matrix->new_from_bioperl($align);
	my $nchar  = $matrix->get_nchar;
	ok( $nchar >= 1140, "$nchar >= 1140 using $tool");
}