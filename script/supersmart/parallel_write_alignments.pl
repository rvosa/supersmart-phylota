#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Matrices::Matrix;
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

# instantiate nodes from infile
my @nodes = $mts->get_nodes_for_table( '-file' => $infile ); ##sequential { return $mts->get_nodes_for_table( '-file' => $infile ) };

# this is sorted from more to less inclusive
my @sorted_clusters = $mts->get_clusters_for_nodes(@nodes); #sequential { return $mts->get_clusters_for_nodes(@nodes) };

# hear we split the sorted clusters into subsets that divide
# the inclusiveness more evenly
my @subset;
my $nworkers = num_workers();
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

my $counter = 0;
my $total = scalar(@clusters);
# make the alignments in parallel mode
my @result = pmap {        
        my $cl = $_;
        my $ti = \%ti;
		$counter = $counter + 1;
		$log->info("processing cluster $counter  / $total");
        my %h = %$cl;             
        # get cluster id
        my $ci = $h{'ci'};
        my @res = ();
        # fetch ALL sequences for the cluster, reduce data set
        my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
        my @ss = $sg->get_sequences_for_cluster_object($cl);
        my @seqs    = $sg->filter_seq_set($sg->get_sequences_for_cluster_object($cl));
        my $single  = $sg->single_cluster($cl);
        my $seed_gi = $single->seed_gi;
        my $mrca    = $single->ti_root->ti;
		$log->info("fetched ".scalar(@seqs)." sequences for cluster $ci");                
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
                                
                my $filename = $workdir . '/' . $seed_gi . '.fa';
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
