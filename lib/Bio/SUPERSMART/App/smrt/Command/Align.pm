package Bio::SUPERSMART::App::smrt::Command::Align;

use strict;
use warnings;

use Bio::Phylo::Matrices::Matrix;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Service::ParallelService 'pthreads';

use Bio::PrimarySeq;
use Bio::AlignIO;

use base 'Bio::SUPERSMART::App::SubCommand';
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
	["outfile|o=s", "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],	
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
    
    # create empty output file 
    open my $outfh, '>', $outfile or die $!;
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
    my @sorted_clusters = $mts->get_clusters_for_nodes(@nodes); 
    my @clusters;	
    my $nworkers = num_workers();	
    
    if (scalar @clusters >= $nworkers ) {
	my @subset;	
	for my $i ( 0 .. $#sorted_clusters ) {
	    my $j = $i % $nworkers;
	    $subset[$j] = [] if not $subset[$j];
				push @{ $subset[$j] }, $sorted_clusters[$i];
	}		
	push @clusters, @{ $subset[$_] } for 0 .. ( $nworkers -1 );
    }
    else {
	@clusters = @sorted_clusters;
    }	
    
    # this is a simple mapping to see whether a taxon is of interest
    my %ti =  map { $_->ti => 1 } @nodes;        
	
    # collect and filter sequences for all clusters
    my %seqs;
    my %seen;
    $log->info("Going to collect sequences for " . scalar(@clusters) . " clusters");
    foreach my $cl (@clusters) {
	# get cluster information
	my $ci = $cl->{'ci'};
	my $type = $cl->{'cl_type'};
	my $single  = $sg->single_cluster($cl);
	my $seed_gi = $single->seed_gi;
	my $mrca    = $single->ti_root->ti;
	
	# make string to identify cluster
	my $clinfo =  $seed_gi . '-' . $mrca . '-' . $ci . '-' . $type;
	my @seqs    = $sg->filter_seq_set($sg->get_sequences_for_cluster_object($cl));
	
	$log->debug("fetched " . scalar(@seqs) . " sequences for cluster id $ci");                
	
	my $t = $seqs[0]->ti;
	
	my @matching = grep { $ti{$_->ti} } @seqs;
	
	# filter out sequences that we have processed before
	@matching = grep { ! $seen{$_->gi} } @matching;
	    
	# let's not keep the ones we can't build alignemnts from
	if ( scalar @matching > 1 ) {
	    foreach my $s (@matching) {
		$seen{$s->gi} = 1;			    
	    }
	    # make sequence objects
	    my @convertedseqs;
	    for my $seq(@matching){						
		my $gi = $seq->gi;
		my $ti  = $sg->find_seq($gi)->ti;
		my $seqobj = Bio::PrimarySeq->new( -display_id => "gi|${gi}|seed_gi|${seed_gi}|taxon|${ti}|mrca|${mrca}" ,
						   -seq => $seq->seq, 
						   -name => $gi,
						   -type => 'dna');
		push @convertedseqs,$seqobj;		    
	    }
	    # store sequences for this cluster
	    $seqs{$clinfo} = \@convertedseqs;
	}	    
    }
    
	# now make the alignments in parallel mode
    my @result = pmap {        		
	my $clinfo = $_;
        
	# alignments could be pre-computed so don't align again if fasta file already exists 
	my $filename = $self->workdir . '/' . $clinfo  . '.fa';	    
	if ( not ( -e  $filename ) ) {	
	    
	    my $s = $seqs{$clinfo};
	    my $aligner = $sg->_make_aligner;
	    $log->info("going to align " . scalar(@{$s}) . " sequences");
	    my $aln = $aligner->align($s);
	    
	    # write alignemnt to fasta file  	                
	    my $out = Bio::AlignIO->new(
		-file => ">$filename",
		-format => 'fasta');
	    $out->write_aln($aln);	                
	} 	
	else {
	    $log->debug("alignment with name $filename already exists, skipping");
	}
	#print alignment file name to output file so it can be saved in output file of script          
	open my $outfh, '>>', $outfile or die $!;
	print $outfh $filename . "\n";
	close $outfh;
	    
    } keys(%seqs);
    
    $log->info("DONE, results written to $outfile");	
    
    return 1;
}

1;
