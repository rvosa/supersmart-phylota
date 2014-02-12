#!/usr/bin/perl
use strict;
use warnings;
use Bio::Phylo::Matrices::Matrix;
use Test::More 'no_plan';

# this tests if we can use the SequenceGetter package
use_ok('Bio::Phylo::PhyLoTA::Service::SequenceGetter');

# this just tests that we can create a new SequenceGetter object
my $sg=new_ok('Bio::Phylo::PhyLoTA::Service::SequenceGetter');

# this fetches the smallest containing cluster around the seed
# sequence with gi 326632174, which is a lemur cytb sequence
my @sequences=$sg->get_smallest_cluster_for_sequence(326632174);

# this computes the most occurring (i.e. median) sequence length,
# which for cytochrome B in primates is 1140 basepairs
my $median = $sg->compute_median_seq_length(@sequences);
ok($median == 1140);

# reduces the number of sequences in the seq by removing duplicate
# sequences and by preferentially retaining sequences of the median
# length, otherwise longer sequences and then otherwise short sequences
my @filtered = $sg->filter_seq_set(@sequences);

# align using muscle
my $align = $sg->align_sequences(@filtered);

# converts a Bio::AlignI object (as returned by muscle) into
# a matrix object that can be written to nexus
my $matrix = Bio::Phylo::Matrices::Matrix->new_from_bioperl($align);
ok($matrix->get_nchar == 1140);

# this part is a brief song and dance to make mesquite happy: a nexus
# file according to mesquite must start with a taxa block and the rows
# in the subsequence matrix must be in the same order as the taxa in the
# taxa block. The make_taxa method generates such a taxa block, but in
# alphabetical order, so we also need to sort the rows in our matrix
# in the same way.
my $taxa = $matrix->make_taxa;
my @sorted = sort { $a->get_name cmp $b->get_name } @{ $matrix->get_entities };
$matrix->clear;
$matrix->insert($_) for @sorted;

# now print out the nexus statements
print "#NEXUS\n",$taxa->to_nexus,$matrix->to_nexus;