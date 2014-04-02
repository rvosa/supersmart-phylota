#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Service::ParallelService 'mpi' => 4; # can be either 'pthreads' or 'mpi';

=head1 NAME

parallel_write_taxa_table.pl - reconciles taxon labels with a taxonomy

=head1 SYNOPSYS

 $ mpi_write_taxa_table.pl --infile=<file> > <outfile>

=head1 DESCRIPTION

Given a file with putative taxon names, attempts to link these to the pipeline's
underlying taxonomy. Writes the results as a tab-separated file to STDOUT.

=cut

# process command line arguments
my $verbosity = WARN;
my @levels = qw[species genus family order class phylum kingdom];
my ( $infile );
GetOptions(
	'infile=s' => \$infile,
	'verbose+' => \$verbosity,
	'level=s'  => \@levels,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => [qw(
		main Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector

	)]
    );

my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;

# read names from file or STDIN, clean line breaks
my @names;
if ( $infile eq '-' ) {
	@names = <STDIN>;
	chomp(@names);
	$log->info("read species names from STDIN");
}
else {
	open my $fh, '<', $infile or die $!;
	@names = <$fh>;
	chomp(@names);
	$log->info("read ".scalar(@names)." species names from $infile");
}

@names = $mts->expand_taxa( @names );

# print table header to stdout in sequential mode
sequential { print join ("\t", 'name', @levels), "\n"; };


# this will take some time to do the taxonomic name resolution in the
# database and with webservices. The below code runs in parallel
my @result = pmap {
        my $name = $_;                
        my @res = ();
        my @nodes = $mts->get_nodes_for_names($name);
        if ( @nodes ) {
                if ( @nodes > 1 ) {
                        $log->warn("found more than one taxon for name $name");
                }
                else {
                        $log->info("found exactly one taxon for name $name");
                }
                
                # for each node, fetch the IDs of all taxonomic levels of interest
                for my $node ( @nodes ) {
                        
                        # create hash of taxonomic levels so that when we walk up the
                        # taxonomy tree we can more easily check to see if we are at a
                        # level of interest
                        my %level = map { $_ => undef } @levels;
                        
                        # traverse up the tree
                        while ( $node ) {
                                my $tn = $node->taxon_name;
                                my $ti = $node->ti;
                                my $rank = $node->rank;
                                if ( exists $level{$rank} ) {
                                        $level{$rank} = $node->get_id;
                                }
                                $node = $node->get_parent;
                        }
                        print join("\t", $name, @level{@levels})."\n";
                        push @res, join("\t", $name, @level{@levels});
                }
        }
        else {
                $log->warn("couldn't resolve name $name");
        }                
        return  @res;
} @names; 
