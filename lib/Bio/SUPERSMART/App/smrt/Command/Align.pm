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
    
    #  If the infile is absent or empty, abort  
    my $file = $opt->infile;
    $self->usage_error("no infile argument given") if not $file;
    $self->usage_error("file $file does not exist") unless (-e $file);
    $self->usage_error("file $file is empty") unless (-s $file);
    
    # Outfile must be empty because we will append to it
    my $outfile = $opt->outfile;
    if ( -e $outfile and -s $outfile ) {
        $self->logger->warn("$outfile is not empty. Will overwrite.");
        open my $outfh, '>', $outfile or die $!;
        close $outfh;       
    }   
}

sub run {
    my ($self, $opt, $args) = @_;
    
    # collect command-line arguments
    my $infile  = $opt->infile;
    my $outfile = $opt->outfile;
    
    # create directory for alignments if specified in argument
    my $dir = $self->outputdir($opt->dirname);  
    
    # instantiate helper objects
    my $log = $self->logger;
    my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
    my $mt  = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
    my $sg  = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
    my $ps  = 'Bio::Phylo::PhyLoTA::Service::ParallelService';
    
    # instantiate taxonomy nodes from infile
    my @nodes = $mts->get_nodes_for_table($mt->parse_taxa_file($infile));     
    $log->info("Found ".scalar(@nodes)." nodes in taxa file $infile");
    
    # fetch all clusters for the focal nodes and organize them into
    # roughly even-sized chunks (they were initially sorted from big
    # to small so the first worker would get all the big ones).
    my @sorted_clusters = $mts->get_clusters_for_nodes(@nodes); 
    my @clusters = $ps->distribute(@sorted_clusters);   
    
    # NOTE: We do the below in sequential mode, which may take some time. This
    # is for the following reason: One sequence (GI) can easily be part of
    # multiple clusters. Since including one sequence into multiple alignments
    # (which would then most likely be merged later anyway ) uses space and
    # resources, we do not include a GI  more than once. To make sure that we
    # skip the GI always in the same alignment, we avoid race conditions in
    # parallel processing. Also see issue #56
    # Collect and filter sequences for all clusters
    my %ti =  map { $_->ti => 1 } @nodes;
    my ( %seqs, %seen, @clinfos );
    $log->info("Going to collect sequences for " . scalar(@clusters) . " clusters");
    for my $cl ( @clusters ) {
    
        # get cluster information
        my $ci      = $cl->{'ci'};
        my $type    = $cl->{'cl_type'};
        my $single  = $sg->single_cluster($cl);
        my $seed_gi = $single->seed_gi;
        my $mrca    = $single->ti_root->ti;
        my $clinfo  = $single->clinfo;
        
        # XXX we should be able to parameterize this more so that we can filter
        # a set down to $max haplotypes per taxon
        my @seqs = $sg->filter_seq_set($sg->get_sequences_for_cluster_object($cl));   
        $log->debug("fetched ".scalar(@seqs)." sequences for cluster $clinfo");
    
        # filter out sequences that we have processed before or that are
        # of uninteresting taxa
        my @matching = grep { $ti{$_->ti} } @seqs;
        @matching    = grep { ! $seen{$_->gi} } @matching;
        
        # let's not keep the ones we can't build alignments from
        if ( scalar @matching > 1 ) {
            $log->info("Will align cluster: $clinfo");
            $seqs{$clinfo} = [
                map {
                    $seen{$_->gi} = 1;          # skip this seq next time                                       
                    $_->to_primary_seq(         # Bio::PrimarySeq
                        'mrca'    => $mrca,
                        'seed_gi' => $seed_gi,
                    );
                } @matching
            ];
            push @clinfos, $clinfo; # to keep optimized distribution for pmap
        }       
    }

    # now make the alignments in parallel mode
    my @result = pmap {             
        my $clinfo = shift;
        
        # align to file, log result
        my $filename = File::Spec->catfile( $dir, $clinfo.'.fa' );      
        $sg->align_to_file( $seqs{$clinfo} => $filename );         
        if ( -s $filename ) {
            $log->info("Have alignment file: $filename");
            open my $outfh, '>>', $outfile or die $!;
            print $outfh $filename . "\n";
            close $outfh;
        }
        
    } @clinfos;
    
    $log->info("DONE, results written to $outfile");    
    
    return 1;
}

sub outputdir {
    my ( $self, $dirname ) = @_;
    my $log = $self->logger;
    
    # create directory for alignments if specified in argument
    my $dir;
    if ( my $alndir = $dirname ) {
        
        # is absolute
        if ( $alndir =~ /^\// ) {
            $dir = $alndir . '/';
        }
        
        # relative
        else {
            $dir = $self->workdir . '/' . $alndir . '/'; 
        }
        $log->info("creating directory $dir");
        mkdir $dir or die $! if not -e $dir;  
    }
    else {
        $dir = $self->workdir . '/'; 
    }
    $log->debug("directory to save alignments in : $dir");
    return $dir;
}
    
1;
