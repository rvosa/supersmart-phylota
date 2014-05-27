#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Parallel::MPI::Simple;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Matrices::Matrix;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use constant HEAD_NODE => 0;
use constant CLUSTER_SUBSET => 1;
use constant ALIGNMENTS_SUBSET => 2;
use constant NODE_MAP => 3;

=head1 NAME

mpi_write_alignments.pl - writes phylogenetically informative clusters for taxa

=head1 SYNOPSYS

 $ mpi_write_alignments.pl --infile=<taxa table> --workdir=<dir> > <outfile>

=head1 DESCRIPTION

Given a table of reconciled taxon names (with NCBI taxon IDs), writes all the 
phylogenetically informative clusters for those taxa as multiple sequence alignments
to the working directory. Prints the produced file names to STDOUT, to be 
re-directed to a file.

=cut

# process command line arguments
my $verbosity = WARN;
my ( $infile, $workdir );
GetOptions(
	'infile=s'  => \$infile,
	'workdir=s' => \$workdir,
	'verbose+'  => \$verbosity,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => [qw(
		main
		Bio::Phylo::PhyLoTA::Service::SequenceGetter
		Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector
	)]
);

# here we start in parallel mode
MPI_Init();
my $rank = MPI_Comm_rank(MPI_COMM_WORLD);

# the following block is executed by the head node
if ( $rank == 0 ) {
        
    # MPI_Comm_size returns the total number of nodes. Because we have one
    # head node this needs to be - 1
    my $nworkers = MPI_Comm_size(MPI_COMM_WORLD) - 1;
    $log->info("we have $nworkers nodes available");
    
    # instantiate helper object
    my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;

    # instantiate nodes from infile
    my @nodes = $mts->get_nodes_for_table( '-file' => $infile );

    # this is sorted from more to less inclusive
    my @sorted_clusters = $mts->get_clusters_for_nodes(@nodes);
    
    # hear we split the sorted clusters into subsets that divide
    # the inclusiveness more evenly
    my @subset;
    for my $i ( 0 .. $#sorted_clusters ) {
		my $j = $i % $nworkers;
		$subset[$j] = [] if not $subset[$j];
		push @{ $subset[$j] }, $sorted_clusters[$i];
    }
    
    # now we flatten the subsets again
    my @clusters;
    push @clusters, @{ $subset[$_] } for 0 .. ( $nworkers - 1 );

    # this is a simple mapping to see whether a taxon is of interest
    my %ti = map { $_->ti => 1 } @nodes;        
    
    # for each worker we dispatch an approximately equal chunk of clusters.
    my $nclusters = int( scalar(@clusters) / $nworkers );
    $log->info("each worker will process about $nclusters clusters");
    
    # this is the starting index of the subset
    my $start = 0;
    
    # iterate over workers to dispatch jobs
    for my $worker ( 1 .. $nworkers ) {
        
        # first send over the map of nodes to keep
		MPI_Send(\%ti,$worker,NODE_MAP,MPI_COMM_WORLD);
        
		# the ending index of the subset is either the starting index + chunk size or the
		# highest index in the array, whichever fits
		my $end  = ( $start + $nclusters ) > $#clusters ? $#clusters : ( $start + $nclusters );
	
		# create subset
		my @subset = @clusters[ $start .. $end ];
		$log->info("dispatching ".scalar(@subset)." clusters (index: $start .. $end) to worker $worker");
		MPI_Send(\@subset,$worker,CLUSTER_SUBSET,MPI_COMM_WORLD);
	
		# increment starting index
		$start += $nclusters + 1;
    }
    
    # iterate over workers to receive and write results
    my $i = 1;
    for my $worker ( 1 .. $nworkers ) {
	
		# get the result
		my $result = MPI_Recv($worker,ALIGNMENTS_SUBSET,MPI_COMM_WORLD);
		$log->info("received ".scalar(@{$result})." results from worker $worker");
	
		# iterate over alignments
		for my $alignment ( @{ $result } ) {
            
            # create out file name
		    my $outfile = $workdir
				. '/'
				. $alignment->{seed_gi}
#				. '.'
#				. $alignment->{gene}
				. '.fa';
            
            # print name to stdout so we can make a list of produced files
            print $outfile, "\n";
            
            # open write handle
		    open my $outfh, '>', $outfile or die $!;
	    
            # iterate over rows in alignment
            for my $row ( @{ $alignment->{matrix} } ) {
                
                # 0 is FASTA header, 1 is aligned sequence data
                print $outfh $row->[0], "\n", $row->[1], "\n";
            }
            
            $i++;
		}
    }    
}

# this block is executed by worker nodes
else {
    
    # first receive the map of nodes to keep
    my $ti = MPI_Recv(0,NODE_MAP,MPI_COMM_WORLD);
    
    # then receive the list of clusters to process
    my $subset = MPI_Recv(0,CLUSTER_SUBSET,MPI_COMM_WORLD);
    
    # iterate over clusters, build result
    my @result;
    for my $cl ( @{ $subset } ) {
        my %h = %$cl;             
        my $ci = $h{'ci'};
        $log->info("worker $rank going to fetch sequences for cluster $ci");
        
        # fetch ALL sequences for the cluster, reduce data set
        my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
        my @seqs    = $sg->filter_seq_set($sg->get_sequences_for_cluster_object($cl));
        my $single  = $sg->single_cluster($cl);
        my $seed_gi = $single->seed_gi;
		my $mrca    = $single->ti_root->ti;
        $log->info("worker $rank fetched ".scalar(@seqs)." sequences");
        
        # keep only the sequences for our taxa
        my @matching = grep { $ti->{$_->ti} } @seqs;
        
        # let's not keep the ones we can't build trees out of
        if ( scalar @matching > 3 ) {
                
            # this runs muscle, so should be on your PATH.
            # this also requires bioperl-live and bioperl-run
            $log->info("worker $rank going to align sequences");
            my $aln = $sg->align_sequences(@matching);
            $log->info("worker $rank done aligning");
            
            # in some cases the alignment tool can crash
            #  here we prevent the whole pipeline to stop in case of a crash
            if ( ! $aln ) {
                    $log->warn("alignment tool crashed! no alignment produced")
            }
            else {
          
                    # convert AlignI to matrix for pretty NEXUS generation
                    my $m = Bio::Phylo::Matrices::Matrix->new_from_bioperl($aln);
                    
                    # iterate over all matrix rows
                    my @matrix;
                    $m->visit(sub{					
                            my $row = shift;
                            
                            # the GI is set as the name by the alignment method
                            my $gi  = $row->get_name;
                            my $ti  = $sg->find_seq($gi)->ti;
                            my $seq = $row->get_char;
                            push @matrix, [ ">gi|${gi}|seed_gi|${seed_gi}|taxon|${ti}|mrca|${mrca}" => $seq ];
                              });
                    
			my $filename = $workdir . '/' . $seed_gi . '.fa';
			open my $fh, '>', $filename or die $!;
			for my $row ( @matrix ) {
				print $fh $row->[0], "\n", $row->[1], "\n";
			}
			
            push @result, {
				'seed_gi' => $seed_gi,
				'matrix'  => \@matrix,
			};
            }
        }        
    }
    
    # return result
    MPI_Send(\@result,HEAD_NODE,ALIGNMENTS_SUBSET,MPI_COMM_WORLD);    
}

# we need to exit cleanly by making this call on all nodes
MPI_Finalize();
