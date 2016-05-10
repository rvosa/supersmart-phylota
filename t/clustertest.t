#!/usr/bin/perl
use strict;
use warnings;
use Bio::SUPERSMART::Service;
use Test::More 'no_plan';

# connect to the database
my $service = Bio::SUPERSMART::Service->new;

# some example GIs:
# matK: 6174789
# ITS1: 18028304
# CO1: 305690971
# rbcL: 149389752
my $seq = $service->find_seq(326632174);

# this is not always going to hold, only for cytochrome B
ok($seq->length == 1140, "found sequence");

# search for all ci/gi intersections that contain focal GI and that are subtrees
my $cigis = $service->search_ci_gi({ 
	'gi'      => $seq->gi,  # sequence ID
	'cl_type' => "subtree", # cluster type
});

# some temporary variables to keep track of cluster size and root taxon 
my ( $clustersize, $biggestcluster, $taxonid ) = ( 0 );

# iterate over all intersections
while( my $c = $cigis->next ) {
    
    # check cluster identifier, root taxon identifier
    ok( $c->clustid, "cluster ID" );
    ok( $c->ti, "taxon ID" );

    # fetch the actual cluster objects
    my $clusters = $service->search_cluster({ 
    	'ci'      => $c->clustid, 
    	'ti_root' => $c->ti,
    });
    
    # iterate over clusters
    while ( my $cluster = $clusters->next ) {
        
        # looking for the largest cluster, i.e. keeping a running tally 
        if ( $cluster->n_ti > $clustersize ){
            ok( $cluster->n_ti > $clustersize, "bigger cluster" );        
            $clustersize = $cluster->n_ti;
            $biggestcluster = $cluster->ci;
            $taxonid = $cluster->ti_root;
        }
    }
}

# now get the biggest cluster
ok( $biggestcluster, "biggest: $biggestcluster" );
ok( $taxonid->ti, "taxon ID of biggest" );

# here's how to fetch the sequences
my $gis = $service->search_ci_gi({ 
	'cl_type' => 'subtree', 
	'ti'      => $taxonid->ti, 
	'clustid' => $biggestcluster
});
ok( $gis, "found sequences" );