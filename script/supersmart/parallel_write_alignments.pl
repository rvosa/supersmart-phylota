#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util 'min';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Matrices::Matrix;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Service::ParallelService 'mpi'; # can be either 'pthreads' or 'mpi';

=head1 NAME

parallel_write_alignments.pl - writes phylogenetically informative clusters for taxa

=head1 SYNOPSYS

 $ parallel_write_alignments.pl --infile=<taxa table> --workdir=<dir> > <outfile>

=head1 DESCRIPTION

Given a table of reconciled taxon names (with NCBI taxon IDs), writes all the 
phylogenetically informative clusters for those taxa as multiple sequence alignments
to the working directory. Prints the produced file names to STDOUT, to be 
re-directed to a file.

=cut

# process command line arguments
my $verbosity = INFO;
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
		Bio::Phylo::PhyLoTA::Service::ParallelService
	)]
);
    
# instantiate helper object
my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
my $mt =  Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
   
# parse the taxa file 
my @taxatable = $mt->parse_taxa_file($infile);

# instantiate nodes from infile
my @nodes = $mts->get_nodes_for_table(@taxatable); 

$log->info("Found " . scalar(@nodes) . " nodes for taxa table");

# this is sorted from more to less inclusive
my @sorted_clusters = $mts->get_clusters_for_nodes(@nodes); 
$log->info("Retrieved clusters for node obects");

# store all the clusters according to their respective seed gi
my %clusters;

sequential{
	foreach my $cl (@sorted_clusters) {
		my $single  = $sg->single_cluster($cl);	
		my $seed_gi = $single->seed_gi;
		#$log->info("Extracted seed gi from cluster");
		push @{$clusters{$seed_gi}}, $single;
	}

	$log->info("Going to reduce clusters");

	# Now reduce the clusters: There can be multiple clusters for one marker (seed gi).
	#  If we have two clusters with the same seed_gi, we want to keep the cluster that
	#  corresponds to the taxon which has higher level, and forget about the other one.
	#  But caution! We do not want a higher level than the highest level in our taxon list!
	my $counter = 0;
	foreach my $seed (keys %clusters) {
		$log->info("Check if clusters for seed gi $seed  can be reduced ( #" . ++$counter . " / " . scalar(keys %clusters) . " )" );
		# collect all possible ranks for clusters for that seed gi
		my @ranks = map { $mts->get_rank_for_taxon( $_->ti_root->ti ) } @{$clusters{$seed}};	
		
		# get all taxonomic ranks ordered from higher to lower levels
		my @all_ranks = $mts->get_taxonomic_ranks;
		# get 'root' rank in taxa table
		my $highest_rank = $mt->get_root_taxon_level( @taxatable );
		# remove everything below the highest taxon
		shift @all_ranks while (not $all_ranks[0] eq $highest_rank);
	
		my @indices;
		foreach my $rank (@ranks) {
			if ($rank ~~ @all_ranks){
				my ($idx) = grep { $all_ranks[$_] eq $rank } 0..$#all_ranks;
				push @indices, $idx;
			}
		}
		
		# Skip the gi if all clusters for this seed gi have ranks lower than the lowest for our taxa
		if (! @indices){
			$log->info("Deleting all (" . scalar(@ranks) . ") cluster(s) for seed gi $seed since they are of higher rank than $highest_rank");
			##$clusters{$seed} =  @{$clusters{$seed}}[0];
			delete $clusters{$seed};
			next;
		}
		# taxonomic ranks are oredered from high to low, we therefore we take the lowest index
		my ($cluster_idx) = grep { $indices[$_] == min(@indices) } 0..$#indices;
		
		# set the chosen cluster as the only one for this seed
		my $cluster = @{$clusters{$seed}}[$cluster_idx];
		$clusters{$seed} = $cluster;
	}
	
	$log->info("Number of clusters : " . scalar(@sorted_clusters));
	@sorted_clusters = values %clusters;
	$log->info("Number of clusters after filtering: " . scalar(@sorted_clusters));
	
};

my @subset;
my $nworkers = num_workers() || 1;
for my $i ( 0 .. $#sorted_clusters ) {
        my $j = $i % $nworkers;
        $subset[$j] = [] if not $subset[$j];
        push @{ $subset[$j] }, $sorted_clusters[$i];
}
    
# now we flatten the subsets again
my @clusters;
push @clusters, @{ $subset[$_] } for 0 .. ( $nworkers - 1 ) ;


# this is a simple mapping to see whether a taxon is of interest
my %ti =  map { $_->ti => 1 } @nodes;        

my $total = scalar(@clusters);
# make the alignments in parallel mode
my @result = pmap {        
        # get 'single' cluster object
        my $single = $_;
        my $ti = \%ti;
		my $ci = $$single{"_column_data"}{"ci"};
		$log->info("processing cluster with id $ci");
		my $cl_type = $$single{"_column_data"}{"cl_type"};
		my $ti_root = $$single{"_column_data"}{"ti_root"};
		
		my %h = ('ci'=>$ci, 'cl_type'=>$cl_type, 'ti_root'=>$ti_root);        
        my @res = ();
        
        # fetch ALL sequences for the cluster, reduce data set
        my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
        my @ss = $sg->get_sequences_for_cluster_object(\%h);
        my @seqs    = $sg->filter_seq_set($sg->get_sequences_for_cluster_object(\%h));        
        my $seed_gi = $single->seed_gi;
        my $mrca    = $single->ti_root->ti;
		$log->info("fetched ".scalar(@seqs)." sequences for cluster id $ci, seed gi: $seed_gi, ti_root : " . $mrca);                
        my @matching = grep { $ti->{$_->ti} } @seqs;
         
        # let's not keep the ones we can't build trees out of
        if ( scalar @matching > 3 ) {
                # this runs muscle or mafft, so should be on your PATH.
                # this also requires bioperl-live and bioperl-run
                $log->info("going to align ".scalar(@matching)." sequences for cluster $ci");
                my $aln = $sg->align_sequences(@matching);
                
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
                                
                my $filename = $workdir . '/' . $seed_gi . '.fa'; #'_' . $mrca . '.fa';
                open my $fh, '>', $filename or die $!;
                for my $row ( @matrix ) {
                        print $fh $row->[0], "\n", $row->[1], "\n";
                }
                close $fh;
                push @res, {
                        'seed_gi' => $seed_gi,
                        'matrix'  => \@matrix,
                };
                $log->info("aligning sequences for cluster $ci finished and written to $filename");
 				
 				#print alignment file name to STDOUT so it can be saved in output file of script          
 				print $filename . "\n";
        } 
        return @res;
} @clusters;
