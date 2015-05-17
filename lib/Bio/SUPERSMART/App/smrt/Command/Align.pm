package Bio::SUPERSMART::App::smrt::Command::Align;

use strict;
use warnings;

use Bio::Phylo::Matrices::Matrix;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Service::ParallelService;

use Bio::PrimarySeq;
use Bio::AlignIO;
use File::Spec;

use Bio::SUPERSMART::App::SubCommand;
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
    my $taxa_default = "species.tsv";
    return (
    ["infile|i=s", "taxa file (tab-seperated value format) as produced by 'smrt taxize'", { arg => "file", default => $taxa_default }],
    ["outfile|o=s", "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],   
    ["dirname|d=s", "write alignments to specified directory name, if not given, alignments are written to working directory", {}]
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

# XXX this method is obscenely long and needs to be broken
# down into more manageable blocks.
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
    my $mt  = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
    my $sg  = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
    
    # instantiate nodes from infile
    my @nodes = $mts->get_nodes_for_table($mt->parse_taxa_file($infile));     
    $log->info("Found " . scalar(@nodes) . " nodes for taxa table");
    
    # fetch all clusters for the focal nodes and organize them into
    # roughly even-sized chunks (they were initially sorted from big
    # to small so the first worker would get all the big ones).
    my @sorted_clusters = $mts->get_clusters_for_nodes(@nodes); 
    my @clusters;   
    my $nworkers = num_workers();
    if ( scalar @sorted_clusters >= $nworkers ) {
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
    
    # collect and filter sequences for all clusters
    # NOTE: We do the below in sequencial mode, which may take some time. This
    # is for the following reason: One sequence (GI) can easily be part of multiple
    # clusters. Since includingh ons sequence into multiple alignments (which would then
    # most likely merged later anyway ) uses space and resources, we do not include a GI 
    # which occured once. To make sure that we skip the GI always in the same alignment,
    # we avoid the race conditions in parallel processing. Also see issue #56

    my %ti =  map { $_->ti => 1 } @nodes;
    my ( %seqs, %seen );
    $log->info("Going to collect sequences for " . scalar(@clusters) . " clusters");
    for my $cl ( @clusters ) {
    
        # get cluster information
        my $ci = $cl->{'ci'};
        my $type = $cl->{'cl_type'};
        my $single  = $sg->single_cluster($cl);
        my $seed_gi = $single->seed_gi;
        my $mrca    = $single->ti_root->ti;
    
        # make string to identify cluster
        my $clinfo =  $seed_gi . '-' . $mrca . '-' . $ci . '-' . $type;
        my @seqs   = $sg->filter_seq_set($sg->get_sequences_for_cluster_object($cl));   
        $log->debug("fetched " . scalar(@seqs) . " sequences for cluster id $ci");                  
        my $t = $seqs[0]->ti;                  
    
        # filter out sequences that we have processed before or that are
        # of uninteresting taxa
        my @matching = grep { $ti{$_->ti} } @seqs;
        @matching    = grep { ! $seen{$_->gi} } @matching;
        
        # let's not keep the ones we can't build alignments from
        if ( scalar @matching > 1 ) {
            for my $s (@matching) {
                $seen{$s->gi} = 1;              
            }
        
            # make sequence objects
            my @convertedseqs;
            for my $seq ( @matching ) {
                my $gi = $seq->gi;
                my $ti  = $sg->find_seq($gi)->ti;
                my $seqobj = Bio::PrimarySeq->new( 
                    '-display_id' => "gi|${gi}|seed_gi|${seed_gi}|taxon|${ti}|mrca|${mrca}",
                    '-seq'        => $seq->seq, 
                    '-name'       => $gi,
                    '-type'       => 'dna'
                );
                push @convertedseqs,$seqobj;            
            }
            # store sequences for this cluster
            $seqs{$clinfo} = \@convertedseqs;
        }       
    }
    
    # create directory for alignments if specified in argument
    my $dir;
    if ( my $alndir = $opt->dirname ) {
	    if ( $alndir =~ /^\// ) {
		    $dir= $alndir . '/';
	    }
	    else {
		    $dir = $self->workdir . '/' . $alndir . '/'; 
	    }
	    $log->info("creating directory $dir");
	    mkdir  $dir or die $! if not -e $dir;  
    }
    else {
	    $dir = $self->workdir . '/'; 
    }
    $log->debug("directory to save alignments in : $dir");

    # now make the alignments in parallel mode
    my @result = pmap {             
        my ($clinfo) = @_;
        
	# create filename
	my $filename = File::Spec->catfile( $dir, $clinfo.'.fa' );
	$log->debug("alignment file name : $filename");

	# alignments could be pre-computed so don't 
	# align again if fasta file already exists  
       if ( not ( -e  $filename ) ) {  
        
            my $s = $seqs{$clinfo};
            my $aligner = $sg->_make_aligner;
            $log->info("going to align " . scalar(@{$s}) . " sequences");
            my $aln = $aligner->align($s);
        
            # write alignment to fasta file                     
            my $out = Bio::AlignIO->new(
                '-file'   => ">$filename",
                '-format' => 'fasta'
            );
            $out->write_aln($aln);                  
        }   
        else {
            $log->debug("alignment with name $filename already exists, skipping");
        }
        
        # print alignment file name to output file 
        # so it can be saved in output file of script          
        open my $outfh, '>>', $outfile or die $!;
        print $outfh $filename . "\n";
        close $outfh;
        
    } keys(%seqs);
    
    $log->info("DONE, results written to $outfile");    
    
    return 1;
}

1;
