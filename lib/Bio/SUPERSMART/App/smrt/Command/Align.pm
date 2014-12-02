package Bio::SUPERSMART::App::smrt::Command::Align;

use strict;
use warnings;

use List::Util 'min';
use Bio::Phylo::Matrices::Matrix;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

use Bio::Seq;
use Bio::AlignIO;

use Parallel::parallel_map;

use base 'Bio::SUPERSMART::App::smrt::SubCommand';

use Bio::SUPERSMART::App::smrt qw(-command);


# ABSTRACT: writes multiple sequence alignments for a list of taxa

=head1 NAME

Align.pm - writes phylogenetically informative clusters for taxa

=head1 SYNOPSYS

smrt align [-h ] [-v ] [-w <dir>] -i <file> [-o <file>] 

=head1 DESCRIPTION

Given an input table of resolved taxa, performs multiple sequence alignment 
of all potentially useful PhyLoTA clusters for the taxa of interest. 
Produces a list of aligned candidate clusters. The alignment method 
can be configured, effective methods include those provided by muscle and mafft.

=cut

sub options {

	my ($self, $opt, $args) = @_;		
	my $outfile_default = "aligned.txt";
	return (
		["infile|i=s", "taxa file (tab-seperated value format) as produced by 'smrt taxize'", { arg => "file", mandatory => 1}],
		["outfile|o=s", "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => "filename"}],	
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	# We only have to check the 'infile' argument. 
	#  If the infile is absent or empty, abort  
	my $file = $opt->infile;
	$self->usage_error("no infile argument given") if not $file;
	$self->usage_error("file $file does not exist") unless (-e $file);
	$self->usage_error("file $file is empty") unless (-s $file);			
}

sub run {
	my ($self, $opt, $args) = @_;

	# collect command-line arguments
	my $infile = $opt->infile;
	my $outfile = $opt->outfile;
	(my $workdir = $opt->workdir) =~ s/\/$//g;

	# create empty output file 
	open my $outfh, '>', $workdir . '/' . $outfile or die $!;
	close $outfh;

	# instantiate helper objects
	my $log = $self->logger;
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $mt =  Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
	   
	# parse the taxa file 
	my @taxatable = $mt->parse_taxa_file($infile);

	# instantiate nodes from infile
	my @nodes = $mts->get_nodes_for_table(@taxatable); 
	
	$log->info("Found " . scalar(@nodes) . " nodes for taxa table");
	
	# this is sorted from more to less inclusive
	my @clusters = $mts->get_clusters_for_nodes(@nodes); 
	
	# this is a simple mapping to see whether a taxon is of interest
	my %ti =  map { $_->ti => 1 } @nodes;        
	
	$log->info("Processing " . scalar(@clusters) . " clusters");
	
	# make the alignments in parallel mode
	my @result = parallel_map {        		
		my $cl = $_;
	   	my $ti = \%ti;
		        
        # returned result: cluster information: seed_gi, mrca, cluster id and -type and finally the alignment
        my @res = ();
        
        # get cluster information
        my $ci = $cl->{'ci'};
        my $type = $cl->{'cl_type'};
        my $single  = $sg->single_cluster($cl);
        my $seed_gi = $single->seed_gi;
        my $mrca    = $single->ti_root->ti;

	   	$log->info("processing cluster with id $ci");
        
        # some alignments can be pre-computed so don't align again if fasta file already exists 
        my $filename = $workdir . '/' . $seed_gi . '-' . $mrca . '-' . $ci . '-' . $type . '.fa';
		
		if ( not ( -e  $filename ) ) {	
	        # fetch ALL sequences for the cluster, reduce data set
	        my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
	        my @ss = $sg->get_sequences_for_cluster_object($cl);
	        my @seqs    = $sg->filter_seq_set($sg->get_sequences_for_cluster_object($cl));
	        
			$log->info("fetched ".scalar(@seqs)." sequences for cluster id $ci");                
	        my @matching = grep { $ti->{$_->ti} } @seqs;
	        
	        # let's not keep the ones we can't build trees out of
	        if ( scalar @matching > 3 ) {
	                
	                # this runs muscle or mafft, so should be on your PATH.
	                # this also requires bioperl-live and bioperl-run
	                                                              
	                $log->info("going to align ".scalar(@matching)." sequences for cluster $ci");
	                my @convertedseqs;
			    	for my $seq(@matching){						
						my $gi = $seq->gi;
						my $ti  = $sg->find_seq($gi)->ti;
						my $seqobj = Bio::Seq->new( -display_id => "gi|${gi}|seed_gi|${seed_gi}|taxon|${ti}|mrca|${mrca}" ,
                             						-seq => $seq->seq, 
                             						-name => $gi,
                             						-type => 'dna');														
						push @convertedseqs,$seqobj;
	     
			    	}
			  
	                my $aligner = $sg->_make_aligner;
    				my $aln = $aligner->align(\@convertedseqs);
	                 
	              	# write alignemnt to fasta file  	                
	                my $out = Bio::AlignIO->new(
	                							-file => ">$filename",
                        		  				-format => 'fasta');
	                $out->write_aln($aln);
	                
	 				#print alignment file name to output file so it can be saved in output file of script          
	 				open my $outfh, '>>', $workdir . '/' . $outfile or die $!;
	 				print $outfh $filename . "\n";
	 				close $outfh;
	        } 
		}
		else {
			$log->debug("alignment with name $filename already exists, skipping");
			# have to report to output file anyways! 
			open my $outfh, '>>', $workdir . '/' . $outfile or die $!;
	 		print $outfh $filename . "\n";
	 		close $outfh;
			
		}
	} @clusters;

	$log->info("DONE, results written to $outfile");
	
	return 1;
}

1;