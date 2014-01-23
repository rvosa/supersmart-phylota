#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

# process command line arguments
my ( $verbosity, $chunksize, $infile ) = ( WARN, 10 );
GetOptions(
    'infile=s'    => \$infile,
    'verbose+'    => \$verbosity,
    'chunksize=i' => \$chunksize,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
    '-level' => $verbosity,
    '-class' => 'main',
);
my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;

# instantiate nodes from infile
my @nodes = $mts->get_nodes_for_table( '-file' => $infile );
my %tip = map { $_->ti => 1 } @nodes;

# first get the most recent common ancestor
$log->info("going to find MRCA for ".scalar(@nodes)." nodes");
my ( $i, %seen, %rank ) = ( 1 );
for my $node ( @nodes ) {
    my @anc = map { $_->ti } @{ $node->get_ancestors };
    
    # ancestors are candidates to be the MRCA if they are
    # an ancestor for every node in the list
    $seen{$_}++ for @anc;
    
    # ancestors are ordered by more to less recent. here
    # we record their rank
    $rank{$anc[$_]} = $_ for 0 .. $#anc;
    $log->info("processed node ".$i++);
}

# these are the nodes that are ancestors of all input nodes
my @candidates = grep { $seen{$_} == scalar(@nodes) } keys %seen;

# this is the shallowest of all shared ancestor nodes, i.e. the MRCA
my ($mrca) = sort { $rank{$a} <=> $rank{$b} } @candidates;
$log->info("MRCA is $mrca");

# now fetch the mrca node object
my $mrca_node = $mts->find_node($mrca);

# do the recursion
recurse($mrca_node);
my ( %terminals, %tips_seen, %parent_of );
sub recurse {
    my $node = shift;
    my $ti = $node->ti;                
    $log->info($node->taxon_name . " ($ti)");
    my @childnodes = @{ $node->get_children };
    recurse($_) for @childnodes;
    my @children = map { $_->ti } @childnodes;
    
    # the node is internal
    if ( @children ) {
        
        # collect all terminal taxa of interest, recursively
        my @all_terminals;
        push @all_terminals, @{ $terminals{$_} } for @children;
        
        # collapse the ones we've already lumped in an earlier internal node
        my @collapsed = map { $parent_of{$_} ? $parent_of{$_} : $_ } @all_terminals;
        my @terminals = sort { $a <=> $b } keys %{{map{$_=>1}@collapsed}};
        $terminals{$ti} = \@terminals;
        $log->info("$ti subtends ".scalar(@terminals)." tips of interest");
        
        # build a unique key for this set of tips
        my $key = join ',', @terminals;
        
        # write the set
        if ( ( scalar(@terminals) > $chunksize and not $tips_seen{$key} ) or $ti == $mrca ) {
            print $ti, "\t", $key, "\n";
            $tips_seen{$key}++;
            $parent_of{$_} = $ti for @terminals;
        }
    }
    
    # the node is terminal
    else {

        # it is a node in the original species list, retain
        if ( $tip{$ti} ) {
            $terminals{$ti} = [ $ti ];
            $log->debug("$ti is one of the focal taxa");
        }

        # not in the original species list, ignore
        else {
            $terminals{$ti} = [];
            $log->debug("ignoring tip $ti");
        }
    }
}