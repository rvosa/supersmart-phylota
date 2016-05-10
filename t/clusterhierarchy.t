#!/usr/bin/perl
use strict;
use warnings;
use Bio::Phylo::Matrices::Matrix;
use Bio::SUPERSMART::Service::SequenceGetter;
use Test::More 'no_plan';

# create sequence getter object
my $sg = Bio::SUPERSMART::Service::SequenceGetter->new;
ok( $sg, "instantiate sequence getter" );

# this fetches the smallest containing cluster around the seed
# sequence with gi 326632174, which is a lemur cytb sequence
my $co = $sg->get_smallest_cluster_object_for_sequence(326632174);
ok( $co, "get small cluster" );

# this fetches the nearest parent cluster
my $parent_co = $sg->get_parent_cluster_object($co);
ok( $parent_co, "get parent cluster" );

# this fetches the nearest children, i.e. we should now have all sibling
# clusters around 326632174
my @child_co = $sg->get_child_cluster_objects($parent_co);
ok( scalar(@child_co), "get sibling clusters" );

for (@child_co) {
    
    # align using muscle
    my $align = $sg->align_sequences($sg->get_sequences_for_cluster_object($_));
    
    # converts a Bio::AlignI object (as returned by muscle) into
    # a matrix object that can be written to nexus
    my $matrix = Bio::Phylo::Matrices::Matrix->new_from_bioperl($align);
    
    ok( $matrix->get_ntax >= 3 );
}