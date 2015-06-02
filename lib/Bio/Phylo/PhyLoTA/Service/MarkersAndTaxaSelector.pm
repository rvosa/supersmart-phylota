# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use strict;
use warnings;
use JSON;
use Moose;
use URI::Escape;
use Data::Dumper;
use LWP::UserAgent;
use Encode 'encode';
use List::MoreUtils qw(firstidx uniq);
use Bio::Phylo::Factory;
use Bio::Phylo::Util::Logger;
use Bio::Phylo::Util::Exceptions 'throw';
use Bio::Phylo::PhyLoTA::DAO;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Service::ParallelService;

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

Accepts a list of name strings, does TNRS on them and returns the equivalent
Node objects from the database.

=cut

sub get_nodes_for_names {
    my ( $self, @names ) = @_;
    
    # Resulting nodes
    my @all_nodes;
    
    # iterate over supplied names
    NAME: for my $name ( @names ) {
        if ( not $name ) {
            $log->warn("Can't search on name '$name', skipping");
            next NAME;
        }
        $log->debug("going to search for name '$name'");
        
        # do we have an exact match?
    # caution: there might be multiple nodes for one taxon name 
        my @nodes = $self->search_node( { taxon_name => $name } )->all;
        
        # no exact match if ->single returns undef (i.e. false)
        if ( scalar @nodes == 0) {
        $log->debug("no exact match for '$name' in local database");
            
            # search the web service
            if ( my $id = $self->_do_tnrs_search($name) ) {
            @nodes = $self->search_node( {  ti => $id } )->all;
            $log->info("Found match $id for $name through TNRS");
        }
        }
        else {
            $log->info("Found exact match(es) for $name in local database");
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
    $log->info("Created " . scalar(@nodes) . " nodes from taxon table");
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
        $log->info("Query completion: $counter/".scalar(@nodes));
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
    eval { 
        my $json_bytes = encode('UTF-8', $result);
        $obj = JSON->new->utf8->decode($json_bytes); 
    
    };
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
                $log->debug("found id $id for name '$name'");
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

=item encode_taxon_name
    
    Encodes a taxon name such that it can be used within SUPERSMART without ambiguity.
    At the moment, we encode by subtituting the following strings/characters:

    '_' is changed to '|', because Bio::Phylo internally converts spaces to '_', which 
    gives trouble with taxon names containing a '_'
    
    ' ' is changed to '_'

=cut

sub encode_taxon_name {
    my ($self, $name) = @_;
    
    # substitute '_' to '|'
    $name =~ s/_/\|/g;
    # substitute ' ' to '_'
    $name =~ s/ /_/g;
    

    return $name;
}

=item decode_taxon_name

    Decodes a taxon name back that was encoded with enceode_taxon_name,
    decoded taxon names can be searched in the database.

=cut

sub decode_taxon_name {
    my ($self, $name) = @_;

    $name =~ s/^'(.*)'$/$1/g;
    $name =~ s/^"(.*)"$/$1/g;                        

    $name =~ s/_/ /g;
    $name =~ s/\|/\_/g;

    return $name;
}

=item make_taxa_table

Given a list of taxon names, returns a taxa table containing 
the taxonomic classifications for each name (if found in database), for all
taxonomic ranks. Returns an array of hashes with keys being ranks (and taxon name), 
values are the taxon IDs at each rank.

=cut

sub make_taxa_table {
    my ($self, @names) = @_;

    # get all possible taxonomic ranks
    my %levels = map{ $_=>1 } $self->get_taxonomic_ranks;
    
    # this will take some time to do the taxonomic name resolution in the
    # database and with webservices. The below code runs in parallel
    my @result = pmap {
        my ($name) = @_;

        # replace consecutive whitespaces with one
        $name =~ s/\s+/ /g;
        my @res   = ();
        my @nodes = $self->get_nodes_for_names($name);
        if (@nodes) {
            if ( @nodes > 1 ) {
                $log->warn("Found more than one taxon for name $name");
            }
            else {
                $log->debug("found exactly one taxon for name $name");
            }
            
            # for each node, fetch the IDs of all taxonomic levels of interest
            for my $node (@nodes) {
                
                # walk up the  taxonomy tree and get tids for levels of interest,
                # store all ids and taxon name in hash
                my %taxinfo;                
                $taxinfo{'name'} = $node->taxon_name;

                # traverse up the tree
                while ( $node ) {
                    my $ti   = $node->ti;
                    my $rank = $node->rank;
                    
                    if ( $levels{$rank} ) {
                        $taxinfo{$rank} = $ti;
                        
                    }
                    $node = $node->get_parent;
                }
                push @res, \%taxinfo;
            }
        }
        else {
            $log->warn("Couldn't resolve name $name");
        }
        return \@res;
    }
    @names;
    
    # flatten result lists (in case there were multiple nodes for one name)
    @result = map {@{$_}} @result;
    $log->info('Created taxa table containing ' . scalar (@result) . ' rows'); 
    
    return @result;
}

=item write_taxa_file

Given a refrence to a taxa table (array of hashes with keys being ranks (and taxon name), values are tids),
writes the file to a tab separated table with headers. Each row corresponds to a taxon, 
first entry is the taxon name followed by the ids for all possible taxonomic ranks, ordered
from low to high. If a taxon ID is not found at a certain level, 'NA' is written into the
table cell.
Omits entries for taxa of higher level than "species".

=cut

sub write_taxa_file {
    my ($self, $outfile, @table) = @_;

    my @levels = reverse ( $self->get_taxonomic_ranks );
    
    open my $out, '>', $outfile or die $!;
    
    # print table header
    print $out join( "\t", 'name', @levels ), "\n";
    
    # only write row for taxon if id for one of the below ranks is given
    my @valid_ranks = ("species", "subspecies", "varietas", "forma");
    
    # loop over table and write each entry, if no id at a taxonomic level exists, write "NA"
    for my $row ( @table ) {
        my %h = %$row;

        # set cell to NA if no taxon ID for rank        
        my @missing_ranks = grep { ! $h{$_}} @levels;
        $h{ $_ } = "NA" for @missing_ranks;
        
        # omit taxa higher than species
        if ( ! grep /[0-9]+/, @h{@valid_ranks} ) {
            $log->info("Omitting higher level taxon " . $h{'name'});
            next;
        }
        print $out join( "\t", $self->encode_taxon_name($h{'name'}),  @h{@levels} ) . "\n";
    }
    close $out;
    $log->info("Wrote taxa table to $outfile");
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

    my @seed_gis;
    
    # parse seed gis from fasta defline
    for my $col (@table) {
        # all entries in one row have the same seed GIs
        my $k = (keys(%$col))[0];
        my $defline = $col->{$k}[0];
        push @seed_gis, $1 if $defline =~ /seed_gi\|([0-9]+)\|/;
    }
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
            if (  my $seqs = $h{$species} ) {
                                my @accessions;
                foreach my $seq ( @$seqs ) {
                    # parse the GI from fasta descripion
                    my ($gi) = $seq =~ /gi\|([0-9]+)\|/;
                    my $seqobj = $self->find_seq($gi);
                    push @accessions,  $seqobj->acc;
                }
                print $outfh join (',', @accessions) . "\t";
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

=item generate_seqids

Returns artificial GI and accession validated to not be in the database already.
These identifiers can be used to insert custom or simulated sequences into the
database. The GI returned is the increment of the highest GI in the database.
The accession is composed of a prefix, an uderscore '_',  and
a random string of length six, consisting of numbers and upper case characters.
Default for prefix is 'SMRT'.

=cut

sub generate_seqids {
    my ($self, $prefix) = @_;   
    my $gi = $self->max_gi + 1; 
    my $acc = $prefix || 'SMRT';    
    my @set= ('0' ..'9', 'A' .. 'Z');
    $acc .= '_' . join ('',map {$set[rand @set]} 1 .. 6);          
    return (acc=>$acc, gi=>$gi);
}


=item enrich_matrix

Given an input matrix, attempts to add up to $config->CLADE_MAX_HAPLOTYPES additional haplotypes to it.

=cut

sub enrich_matrix {
    my ($self, $matrix) = @_;
    my $log = $self->logger;
    my $max = $self->config->CLADE_MAX_HAPLOTYPES;
    my %taxa;
    my %lookup;
    $matrix->visit( sub {

        # cache the existing row
        my $row   = shift;
        my $gi    = $row->get_meta_object("smrt:gi");
        my $ti    = $row->get_taxon->get_meta_object("smrt:tid");
        my $mrca  = $row->get_meta_object("smrt:mrca");
        my $seed  = $row->get_meta_object("smrt:seed_gi");
        my $taxon = $row->get_taxon;
        $matrix->delete($row);
        $taxon->unset_datum($row);
        $taxa{$ti} = {} if not $taxa{$ti};
        $taxa{$ti}->{$gi} = $self->find_seq($gi);
        $lookup{$gi} = { "smrt:tid" => $ti, "smrt:mrca" => $mrca, "smrt:seed_gi" => $seed };

        # find its cluster(s)
        $log->info("Going to search subtree clusters with seed_gi $seed, mrca $mrca");
        my $results = $self->search_cluster({
            'seed_gi' => $seed,
            'ti_root' => $mrca,
            'cl_type' => 'subtree',
        });

        # find other candidates: same species, diff gi
        my ( @candidates );
        while ( my $res = $results->next ) {
            my $ci = $res->ci;
            $log->info("Going to search other haplotypes for species $ti in focal cluster");
            my $joins = $self->search_ci_gi({
                'clustid'  => $ci,
                'ti'       => $mrca,
                'cl_type'  => 'subtree',
                'ti_of_gi' => $ti,
            });
            my ( @c, $seen );
            while ( my $join = $joins->next ) {
                my $cgi = $join->gi;
                if ( $cgi != $gi ) {
                    push @c, $cgi;
                    $log->info("Found candidate: $cgi");
                }
                $seen++ if $cgi == $gi;
            }
            push @candidates, @c if $seen;
        }

        # reduce to MAX_HAPLOTYPES sorted by length
        $log->info("Found ".scalar(@candidates)." sibling haplotypes for $gi");
        my @sorted = sort { $b->length <=> $a->length } map { $self->find_seq($_) } @candidates;
        my %seq = ( $taxa{$ti}->{$gi}->seq => 1 );
        CANDIDATE: while ( $max > scalar keys %{ $taxa{$ti} } ) {
            my $c = shift @candidates;
            if ( $c ) {
                $c = $self->find_seq($c);
                my $seq = $c->seq;
                if ( not $seq{$seq}++ ) {
                    $taxa{$ti}->{$c->gi} = $c;
                    $lookup{$c->gi} = $lookup{$gi};
                }
            }
            last CANDIDATE unless @candidates;
        }
    });

    # now align
    my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
    my @seqs;
    for my $ti ( keys %taxa ) {
        push @seqs, values %{ $taxa{$ti} };
    }
    $log->info("Going to align ".scalar(@seqs)." sequences");
    my $aln = $sg->align_sequences(@seqs);

    # and add to the matrix
    my $taxa = $matrix->get_taxa;
    my $to   = $matrix->get_type_object;
    my $fac  = Bio::Phylo::Factory->new;
    my %counter;
    for my $row ( $aln->each_seq ) {
        my $gi     = $row->id;
        my $seq    = $row->seq;
        my $lookup = $lookup{$gi};
        my $ti     = $lookup->{"smrt:tid"};
        my $suffix = ++$counter{$ti};
        my $taxon  = $taxa->get_by_name( $ti );
        my $datum  = $fac->create_datum(
            '-type_object' => $to,
            '-taxon'       => $taxon,
            '-name'        => $ti . '_' . $suffix,
            '-char'        => $seq,
        );
        $datum->set_meta_object( "smrt:gi" => $gi );
        for my $predicate ( keys %$lookup ) {
            $datum->set_meta_object( $predicate => $lookup->{$predicate} );
        }
        $matrix->insert($datum);
    }
    return $matrix;
}

=item degap_matrix

Given an input matrix, splices out all columns with only gaps

=cut


sub degap_matrix {
    my ( $self, $matrix ) = @_;
    
    my ( %gaps, $ntax );
    $matrix->visit(sub{
        my $row = shift;
        
        # keep track of which columns might be all gaps
        my @char = $row->get_char;
        for my $i ( 0 .. $#char ) {
            my $c = $char[$i];
            if ( $c eq '-' or $c eq '?' or $c eq 'N' ) {
                $gaps{$i}++;                            
            }
        }
        $ntax++;
    });                 
    my $max = $matrix->get_nchar - 1;
    my @indices = grep { not defined $gaps{$_} or $gaps{$_} < $ntax } 0 .. $max;
    if ( @indices < $matrix->get_nchar ) {
        $matrix->visit(sub{
            my $row = shift;
            my @char = $row->get_char;
            my @ungapped = @char[@indices];
            $row->set_char(@ungapped);
        });
        return $matrix;
    }
    else {
        return undef;
    }
}

=begin comment

Private method to retrieve the contents of a URL

=end comment

=cut

# fetch data from a URL
sub _fetch_url {
    my ( $url ) = @_;
    $log->debug("going to fetch $url");
    
    # instantiate user agent
    my $ua = LWP::UserAgent->new;
    $ua->timeout($timeout);
    $log->debug("instantiated user agent with timeout $timeout");
    
    # do the request on LWP::UserAgent $ua
    my $response = $ua->get($url);
    
    # had a 200 OK
    if ( $response->is_success or $response->code == 302 ) {
        $log->debug($response->status_line);
        my $content = $response->decoded_content;
        $log->debug($content);
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
