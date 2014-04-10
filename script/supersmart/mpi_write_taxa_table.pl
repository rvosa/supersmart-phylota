#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::MPI::Simple;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use constant DIRTY_NAMES_SUBSET => 1;
use constant CLEAN_NAMES_SUBSET => 2;
use constant HEAD_NODE => 0;

=head1 NAME

mpi_write_taxa_table.pl - reconciles taxon labels with a taxonomy

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
		main 
	)]
);

# here we start in parallel mode
MPI_Init();
my $rank = MPI_Comm_rank(MPI_COMM_WORLD);

my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;

# this is executed by the boss node
if ( $rank == 0 ) {

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

    # MPI_Comm_size returns the total number of nodes. Because we have one
    # head node this needs to be - 1
    my $nworkers = MPI_Comm_size(MPI_COMM_WORLD) - 1;
    $log->info("we have $nworkers nodes available");
    
    # for each worker we dispatch an approximately equal chunk of names.
    my $nnames = int( scalar(@names) / $nworkers );
    $log->info("each worker will process about $nnames names");
    
    # this is the starting index of the subset
    my $start = 0;
    
    # iterate over workers to dispatch jobs
    for my $worker ( 1 .. $nworkers ) {
	
	# the ending index of the subset is either the starting index + chunk size or the
	# highest index in the array, whichever fits
	my $end  = ( $start + $nnames ) > $#names ? $#names : ( $start + $nnames );
	
	# create subset
	my @subset = @names[ $start .. $end ];
	$log->info("dispatching ".scalar(@subset)." names (index: $start .. $end) to worker $worker");
	MPI_Send(\@subset,$worker,DIRTY_NAMES_SUBSET,MPI_COMM_WORLD);
	
	# increment starting index
	$start += $nnames + 1;
    }
    
    # print table header
    print join("\t", 'name', @levels), "\n";
    
    # iterate over workers to receive results
    for my $worker ( 1 .. $nworkers ) {
	
	# get the result
	my $result = MPI_Recv($worker,CLEAN_NAMES_SUBSET,MPI_COMM_WORLD);
	$log->info("received ".scalar(@{$result})." results from worker $worker");
	
	# print result to STDOUT
	for my $row ( @{ $result } ) {
	    print join( "\t", @{ $row } ), "\n";
	}
    }
}

# this is executed by worker nodes
else {
    my $subset = MPI_Recv(0,DIRTY_NAMES_SUBSET,MPI_COMM_WORLD);
    my @names = @{ $subset };

    # this will take some time to do the taxonomic name resolution in the
    # database and with webservices
    my @result;
    for my $name ( @names ) {
	my @nodes = $mts->get_nodes_for_names($name);
	if ( @nodes ) {
	    if ( @nodes > 1 ) {
		$log->warn("worker $rank found more than one taxon for name $name");
	    }
	    else {
		$log->info("worker $rank found exactly one taxon for name $name");
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
		push @result, [ $name, @level{@levels} ];
	    }
	}
	else {
	    $log->warn("worker $rank couldn't resolve name $name");
	}
    }
    
    # return result
    MPI_Send(\@result,HEAD_NODE,CLEAN_NAMES_SUBSET,MPI_COMM_WORLD);
}

# we need to exit cleanly by making this call on all nodes
MPI_Finalize();
