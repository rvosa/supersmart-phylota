# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use strict;
use warnings;
use JSON;
use Moose;
use URI::Escape;
use Data::Dumper;
use LWP::UserAgent;
use Bio::Phylo::Factory;
use Bio::Phylo::Util::Logger;
use Bio::Phylo::Util::Exceptions 'throw';
use Bio::Phylo::PhyLoTA::DAO;

extends 'Bio::Phylo::PhyLoTA::Service';

=head1 NAME

Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector - Markers and Taxa Selector

=head1 DESCRIPTION

Selects optimal set of taxa and markers based on #markers/sp, coverage on matrix
(total missing data). User can change threshold.

=cut

# URL for the taxonomic name resolution service
my $TNRS_BASE     = 'http://taxosaurus.org/';
my $TNRS_URL      = $TNRS_BASE . 'submit';
my $TNRS_RETRIEVE = $TNRS_BASE . 'retrieve/';

# defaults for web service
my $timeout = 120;
my $wait    = 5;

# this is used to create other objects, i.e. projects and taxa
my $fac = Bio::Phylo::Factory->new;

# this is used for logging messages
my $log = Bio::Phylo::Util::Logger->new;

# these are all ranks in NCBI taxonomy which are considered in SUPERSMART. 
my @taxonomic_ranks = ('superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum',
                               'subphylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder',
                               'order', 'suborder', 'infraorder', 'parvorder', 'superfamily', 'family',
                               'subfamily', 'tribe', 'subtribe', 'genus', 'subgenus', 'species group',
                               'species subgroup', 'species', 'subspecies','varietas', 'forma');


=over

=item get_nodes_for_names

Accepts a list of name string, does TNRS on them and returns the equivalent
Node objects from the database.

=cut

sub get_nodes_for_names {
    my ( $self, @names ) = @_;
    
    # Resulting nodes
    my @all_nodes;
    
    # iterate over supplied names
    NAME: for my $name ( @names ) {
        if ( not $name ) {
            $log->warn("can't search on name '$name', skipping");
            next NAME;
        }
        $log->info("going to search for name '$name'");
        
        # do we have an exact match?
		# caution: there might be multiple nodes for one taxon name 
    	my @nodes = $self->search_node( { taxon_name => $name } )->all;
        
        # no exact match if ->single returns undef (i.e. false)
        if ( scalar @nodes == 0) {
	    $log->info("no exact match for '$name' in local database");
            
            # search the web service
            if ( my $id = $self->_do_tnrs_search($name) ) {
	    		@nodes = $self->search_node( {  ti => $id } )->all;
	    	$log->info("found match $id for $name through TNRS");
	    }
            else {
                $log->warn("couldn't find name $name anywhere!");
            }
        }
        else {
            $log->info("found  exact match(es) for $name in local database");
        }
        
        # store result
        push @all_nodes, @nodes if scalar @nodes > 0;        
    }
    # return results
    return @all_nodes;
}

=item get_nodes_for_table

Accepts an array of taxa table entries as produced by parse_taxa_file in class Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa 
and returns instantiated node objects

=cut

sub get_nodes_for_table {
    my ( $self, @taxatable ) = @_;
   	my @nodes = ();
    
    # iterate over lines
    foreach (@taxatable) {
		
		my %entry = %{$_};
    	
    	# get the entries of the row sorted from lower to higher taxa
        my @fields =  @entry{reverse $self->get_taxonomic_ranks};
		
		# get the id of the lowest taxon that is mapped to the entry and skip e.g. NA values
		my $id;
		foreach my $f ( @fields ) {
			if ( $f ) {
				if ( $f =~ /^\d+$/ ) {
					$id = $f;
					last;
				} 			
			}
		}	
		# omit e.g. header, rows with only NA...
		next if not $id;
		
		my $node = $self->find_node($id);
            
        # this *could* fail because perhaps the IDs have been found by TNRS
        # without being present in our local database
        if ( $node ) {
        	push @nodes, $node;
        }
        else {
            $log->warn("couldn't instantiate node $id");
        }
    }
    $log->info("created " . scalar(@nodes) . " nodes from taxon table");
    return @nodes;
}

=item get_clusters_for_nodes

Given a list of Node objects, returns the clusters that contain these Nodes in
order of decreasing inclusiveness.

Note: if this query runs really, really slowly it is because your version of the
database does not yet have an index on ci_gi.ti_of_gi

To fix this:

mysql> ALTER TABLE ci_gi ADD INDEX(ti_of_gi);

=cut

sub get_clusters_for_nodes {
    my ( $self, @nodes ) = @_;
    my $c = {};
    my @clusters;
    my $counter = 1;
    
    # iterate over nodes
    for my $node ( @nodes ) {
        $log->info("query completion: $counter/".scalar(@nodes));
        $log->debug("finding clusters for ".$node->ti);
        
        # find ci_gi intersection records for this node's ti
        my @cigis = $self->search_ci_gi({ ti_of_gi => $node->ti })->all;
        
        # iterate over matches
        for my $cigi ( @cigis ) {
            
            # store these for re-use: their combination is the composite
            # key to fetch clusters from the clusters table
            my $ti      = $cigi->ti; # this is not the input taxon id, but cluster.ti_root
            my $clustid = $cigi->clustid; # this is the same as cluster.ci
            my $cl_type = $cigi->cl_type;
            $log->debug("ti => $ti, clustid => $clustid, cl_type => $cl_type");
            
            # build path
            $c->{$ti} = {} unless $c->{$ti};
            $c->{$ti}->{$clustid} = {} unless $c->{$ti}->{$clustid};
            
            push @clusters, {
                'ti'      => $ti,
                'clustid' => $clustid,
                'cl_type' => $cl_type,
            } unless $c->{$ti}->{$clustid}->{$cl_type}++;
        }
        $counter++;
    }
    
    # schwartzian transform
    return
    map {{
        'ti_root' => $_->{ti},
        'ci'      => $_->{clustid},
        'cl_type' => $_->{cl_type},
    }}
    sort {
        $b->{cover} <=> $a->{cover}
    }
    map {{
        'ti'      => $_->{ti},
        'clustid' => $_->{clustid},
        'cl_type' => $_->{cl_type},
        'cover'   => $c->{ $_->{ti} }->{ $_->{clustid} }->{ $_->{cl_type} }
    }} @clusters;
}

=item get_tree_for_nodes

Creates the Bio::Phylo tree that spans these nodes

=cut

sub get_tree_for_nodes {
    my ( $self, @nodes ) = @_;
    
    # create new tree
    my $tree = $fac->create_tree;
    $log->debug("created new tree object");
    
    # build tree structure in hashes
    my ( %children, %by_id );
    for my $node ( @nodes ) {
        
        # for each node traverse up to the root
        while( my $parent = $node->get_parent ) {
            my $pid = $parent->get_id;
            my $nid = $node->get_id;
            
            $children{$pid} = {} unless $children{$pid};
            $children{$pid}->{$nid} = 1;
            
            # we will visit the same node via multiple, distinct paths
            if ( not $by_id{$nid} ) {
                #$log->debug("creating node with guid $nid and name ".$node->taxon_name);
                $by_id{$nid} = $fac->create_node( '-guid' => $nid, '-name' => $node->taxon_name );
            }
            $node = $parent;
        }
        
        # this happens the first time, for the root
        my $nid = $node->get_id;
        if ( not $by_id{$nid} ) {
            #$log->debug("creating node with guid $nid and name ".$node->taxon_name);
            $by_id{$nid} = $fac->create_node( '-guid' => $nid, '-name' => $node->taxon_name );
        }        
    }
    
    # copy the tree structure
    for my $nid ( keys %by_id ) {
        my $bio_phylo_node = $by_id{$nid};
        #$log->debug($bio_phylo_node->to_string);
        if ( $children{$nid} ) {
            for my $child_id ( keys %{ $children{$nid} } ) {
                $bio_phylo_node->set_child($by_id{$child_id});
                $by_id{$child_id}->set_parent($bio_phylo_node);
            }
        }
        $tree->insert($bio_phylo_node);
    }
    
    return $tree;
}


=item expand_taxa

Returns a list of taxon names for all descendent leaf taxa of a given root taxon.
The 'lowest taxonomic rank' to be considered can be specified as 
an argument to this function, e.g. from the config file. 
If leaves are found that are lower than the 'lowest taxonomic rank', their parents are
considered a leaf in the tree (if themselfes are at least at the level 
of 'lowest taxonomic rank').   

=cut

sub expand_taxa {
        my ($self, $names, $lowest_rank) = @_;                
        my @root_taxa = @{$names};
    
        # keep ranks that are higher or equal to highest rank specified as argument
        # do not expand taxa if lowest rank not given
        if ( ! $lowest_rank ) {
                $log->info("taxonomic rank for taxa expansion not set, skipping");
                return ( @root_taxa );
        }
        
        $log->info("expanding taxa names to taxonomic level down to '$lowest_rank'");

        my @result = ();

        my ($index) = grep { uc $taxonomic_ranks[$_] eq uc $lowest_rank } 0..$#taxonomic_ranks;
        if (! $index){
                $log->warn("invalid lowest taxon rank '$lowest_rank', considering all ranks for taxa expansion");
        }
        # make subset of ranks that are considered, add 'no rank'
        my @valid_ranks = @taxonomic_ranks[ 0 .. $index ];
        push @valid_ranks, 'no rank';
                
        my @nodes = $self->get_nodes_for_names( @root_taxa );

        my @queue = ( @nodes );

        # traverse the tree upwards and keep the leaves that have at least
        #  the taxnomic rank specified by parameter $lowest_rank
        while ( @queue ) {
                my $node = shift @queue;
                $log->debug("Processing node " . $node->get_name);
                my @children = @{ $node->get_children };

                # if there is at least one child which has not reached the 
                # specified lowest taxonomic rank yet, we keep traversing 
                # the tree
                my @child_ranks = map { $_->rank } @children;
                foreach  my $rank ( @child_ranks ) {
                        if ( grep { $_ eq $rank } @valid_ranks ) {
                				push @queue, @children;                                
                                last;
                        }
                }                
                # check if current node has valid rank, if yes it 
                # is in the result list, except if it is the root taxon that
                #  is expanded!
				my $name = $node->get_name;
                if ( grep { $_ eq $node->rank } @valid_ranks) {  #and not  grep (/^$name$/, @root_taxa)  ) {
                                my $tipname = $name;
                                push @result, $tipname;                        
                }
        }                 
        return @result;
}

=item get_rank_for_taxon

Returns a string with the taxonomic rank for the given taxon id

=cut

sub get_rank_for_taxon {
	my ($self, $tid) = @_;
	return $self->find_node($tid)->rank;
}

=item get_taxonomic_ranks

Getter for all taxonomic ranks which are considered

=cut

sub get_taxonomic_ranks {
	my $self = shift;
	return @taxonomic_ranks;
}

=item get_outgroup_taxa

Given a list of taxon IDs and a tree, returns a list of taxon ids that are an outgroup 
with respect to the given input list. The outgroup is selected by calculating the MRCA of 
of all input species, then select a sister node of the MRCA in the tree and return all its
terminal nodes

=cut

sub get_outgroup_taxa {
	my ( $self, $tree, $ingroup ) = @_;
	
	# get node object for taxon ids (or names) 
	my @ids = @{$ingroup};
	my @nodes;
	$tree->visit(
		sub{ 
			my $node = shift;
			my $name = $node->get_name;
			if (  grep ( $_ eq $name, @ids )) {
				push @nodes, $node; 
			}
		}
	);
	
	my $mrca = $tree->get_mrca(\@nodes);
		
	my @terminals;
	while ( (scalar (@terminals)) == 0 and (! $mrca->is_root) ){
		my $sister = $mrca->get_next_sister || $mrca->get_previous_sister;
		if ($sister) {
			@terminals = @{ $sister->get_terminals };
		} else {
			$mrca = $mrca->get_parent;
		}
	}
	
	if ($mrca->is_root){
		$log->warn("cannot determine outgroup taxa since the mrca of the ingroup is the root of the tree!");	
		return ();	
	}
	
	return ( map {$_->get_name} @terminals );
}


=item taxa_are_disjoint

Computes whether two array references of taxon IDs are non-overlapping (i.e.
disjoint).

=cut

sub taxa_are_disjoint {
    my ( $self, $set1, $set2 ) = @_;
    my %set1 = map { $_ => 1 } map { ref $_ ? $_->ti : $_ } @{ $set1 };
    my %set2 = map { $_ => 1 } map { ref $_ ? $_->ti : $_ } @{ $set2 };
    
    # check if any taxon from set1 occurs in set2
    for my $t1 ( keys %set1 ) {
        if ( $set2{$t1} ) {
            return 0;
        }
    }
    
    # check if any taxon from set2 occurs in set1
    for my $t2 ( keys %set2 ) {
        if ( $set1{$t2} ) {
            return 0;
        }
    }
    
    # the sets are disjoint, so return true
    return 1;
}

=begin comment

Private method for querying the TNRS web service

=end comment

=cut

sub _do_tnrs_search {
    my ( $self, $name ) = @_;
    
    # do the request
    my $result = _fetch_url( $TNRS_URL . '?query=' . uri_escape( $name ) );
    return if not $result;
    $log->debug("raw result: $result");
    
    # sometimes there are UTF-8 encoding errors, which the JSON parser
    # chokes on. the best we can do is try to trap these and give up.
    my $obj;
    eval { $obj = decode_json($result); };
    if ( $@ ) {
    	$log->warn( "error decoding JSON $result: $@" );
    	return;
    }    
    
    $log->debug("parsed response: ".Dumper($obj));
    
    # start polling
    while(1) {
        
        # we have a final response
        if ( $obj->{'names'} ) {
            if ( my $id = $self->_process_matches($obj->{'names'}) ) {
                $log->info("found id $id for name '$name'");
                return $id; # done
            }
            else {
                return;
            }
        }
        
        # do another cycle
        sleep $wait;
        
        # try to reconstruct the retrieve URL. This seems to be variable.
        my $url;
        if ( $obj->{'uri'} ) {
            $url = $obj->{'uri'};
        }
        elsif ( $obj->{'message'} =~ /Job (\S+) is still being processed./ ) {
            $url = $TNRS_RETRIEVE . $1;
        }
        else {
            die "Don't know how to continue";
        }
        
        # try to fetch the retrieve URL. If this fails we can only return undef.
        if ( my $result = _fetch_url($url) ) {
            $obj = decode_json($result);
        }
        else {
            $log->error("No result for $name");
            return;
        }
    }    
}

=begin comment

Private method for processing the JSON returned by TNRS

=end comment

=cut

sub _process_matches {
    my ( $self, $names ) = @_;
    for my $name ( @{ $names } ) {
        for my $match ( @{ $name->{'matches'} } ) {
            if ( $match->{'sourceId'} eq 'NCBI' and $match->{'uri'} =~ /(\d+)$/ ) {
                my $id = $1;
                return $id;
            }
        }
    }
    return;
}


=item write_marker_summary

Writes a table containing all species as rows and all chosen markers  as columns, 
reports the genbank accession of a specific marker for a species.

=cut
sub write_marker_summary {
	my ( $self, $file, $tab, $specs ) = @_;
	my @table = @$tab;
	my @all_species = @$specs;
	
	# remove empty columns (for markers that are never included)
	@table = grep {keys %{$_}} @table;
	
	my @seed_gis = map { (values(%{$_}))[0] =~ /seed_gi\|([0-9]+)\|/ } @table;
	
	open my $outfh, '>', $file or die $!;
	
	# print table header
	print $outfh "taxon\t";
	print $outfh "marker" . $_ ."\t" for 1..$#seed_gis+1;
	print $outfh "\n";
	foreach my $species ( @all_species ) {
		my $name = $self->find_node($species)->taxon_name;
		print $outfh $name . "\t";
		foreach my $column ( @table ) {
			my %h = %{$column};
			if (  my $seq = $h{$species} ) {
				# parse the GI from fasta descripion
				my ($gi) = $seq =~ /gi\|([0-9]+)\|/;
				my $seqobj = $self->find_seq($gi);
				print $outfh $seqobj->acc . "\t";
			}
			else {
				print $outfh "\t";
			}
		}
		print $outfh "\n";
	}
	print $outfh "\n";
	
	# print information about markers at the bottom of the table
	foreach my $i ( 1..$#seed_gis+1 ) {
		my $seqobj = $self->find_seq( $seed_gis[$i-1] );
		print $outfh "# marker$i cluster seed: " . $seqobj->acc . ", description: " . $seqobj->def . "\n"; 
	}
	close $outfh;
}



=begin comment

Private method to retrieve the contents of a URL

=end comment

=cut

# fetch data from a URL
sub _fetch_url {
	my ( $url ) = @_;
	$log->info("going to fetch $url");
	
	# instantiate user agent
	my $ua = LWP::UserAgent->new;
	$ua->timeout($timeout);
	$log->info("instantiated user agent with timeout $timeout");
	
	# do the request on LWP::UserAgent $ua
	my $response = $ua->get($url);
	
	# had a 200 OK
	if ( $response->is_success or $response->code == 302 ) {
		$log->info($response->status_line);
		my $content = $response->decoded_content;
		$log->info($content);
		return $content;
	}
	else {
		$log->error($url . ' - ' . $response->status_line);
		return undef;
	}	
}

=back 

=cut

1;
