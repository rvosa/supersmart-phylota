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
use Bio::Phylo::PhyLoTA::Service::ParallelService 'pthreads'; # can be either 'pthreads' or 'mpi';

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
		Bio::Phylo::PhyLoTA::Service::ParallelService
	)]
);

    
# instantiate helper object
my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;

# instantiate nodes from infile
my @nodes = $mts->get_nodes_for_table( '-file' => $infile ); ##sequential { return $mts->get_nodes_for_table( '-file' => $infile ) };

# this is sorted from more to less inclusive
my @sorted_clusters = $mts->get_clusters_for_nodes(@nodes); #sequential { return $mts->get_clusters_for_nodes(@nodes) };


# this is a simple mapping to see whether a taxon is of interest
my @ti = map { $_->ti } @nodes;

# make the alignments in parallel mode
my @result = pmap {        
        my $cl = $_;
        my @res = ();
        # fetch ALL sequences for the cluster, reduce data set
        my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
        my @ss = $sg->get_sequences_for_cluster_object($cl);
        my @seqs    = $sg->filter_seq_set($sg->get_sequences_for_cluster_object($cl));
        
        my $single  = $sg->single_cluster($cl);
        my $seed_gi = $single->seed_gi;
        my $mrca    = $single->ti_root->ti;
                
        my @matching = grep { my $s = $_; grep { $_ == $s->ti } @ti  } @seqs;
      
        # let's not keep the ones we can't build trees out of
        if ( scalar @matching > 3 ) {
                
                # this runs muscle, so should be on your PATH.
                # this also requires bioperl-live and bioperl-run
                #$log->info("going to align ".scalar(@matching)." sequences");
                my $aln = $sg->align_sequences(@matching);
                #$log->info("done aligning");
            
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
                
                push @res, {
                        'seed_gi' => $seed_gi,
                        'matrix'  => \@matrix,
                };
        }
        return @res;
} @sorted_clusters;


# process results from parallel run
sequential {    
        # iterate over alignments
        for my $alignment ( @result ) {            
                # create out file name
                my $outfile = $workdir
                    . '/'
                    . $alignment->{seed_gi}
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
        }    
};

