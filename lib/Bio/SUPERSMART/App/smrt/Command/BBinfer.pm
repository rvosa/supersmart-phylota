package Bio::SUPERSMART::App::smrt::Command::BBinfer;
use strict;
use warnings;
use Bio::Phylo::IO qw(parse parse_tree);
use Bio::Phylo::Factory;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Service::InferenceService;
use Bio::Phylo::Util::Exceptions 'throw';
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: inference of a genus-level backbone tree

=head1 NAME

BBinfer.pm - inference of genus-level backbone tree

=head1 SYNOPSIS

smrt bbinfer [-h ] [-v ] [-w <dir>] -s <file> -t <file> [-i <tool>] [-o <file>] 

=head1 DESCRIPTION

Given an input supermatrix (in interleaved PHYLIP format), infers a backbone tree or set
of trees. In the case of maximum likelihood inference, the set of trees are bootstrapped
replicates (this applies to the tools ExaML and RAxML, currently). For ExaBayes, this
is sample of the posterior distribution of trees.

The tree(s) resulting from this command are written to file. ExaML and ExaBayes produce 
many intermediate checkpoint files, for which a directory location needs to be specified. 

=cut

# define additional command line arguments
sub options {
    my ($self, $opt, $args) = @_;       
    my $outfile_default = "backbone.dnd";
    my $tool_default = "examl";
    my $boot_default = 1;
    return (
        ["supermatrix|s=s", "matrix of concatenated multiple sequece alignments as produced by 'smrt bbmerge'", { arg => "file", mandatory => 1}],  
        ["starttree|t=s", "starting tree for ExaML tree inference. If not given, a random starting tree is generated", { arg => "file", mandatory => 1}],
        ["inferencetool|i=s", "software tool for backbone inference (RAxML, ExaML or ExaBayes), defaults to $tool_default", {default => $tool_default, arg => "tool"}],         
        ["bootstrap|b=i", "number of bootstrap replicates. Will add the support values to the backbone tree. Not applicable to Bayesian methods.", { default => $boot_default }],
        ["ids|n", "Return tree with NCBI identifiers instead of taxon names", {}],      
        ["outfile|o=s", "name of the output tree file (in newick format), defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],
    );  
}

# validate command line arguments
sub validate {
    my ($self, $opt, $args) = @_;

    my $sm = $opt->supermatrix;
    my $st = $opt->starttree;
    my $tool = $opt->inferencetool;
    
    $self->usage_error("need supermatrix") if not $sm;
    $self->usage_error("file $sm does not exist") unless -e $sm;
    $self->usage_error("file $sm is empty") unless -s $sm;
}

# run the analysis
sub run {
    my ( $self, $opt, $args ) = @_;     
        
    # collect command-line arguments
    my $supermatrix = $opt->supermatrix;
    my $starttree   = $opt->starttree;
    my $bootstrap   = $opt->bootstrap;
                
    # instantiate and configure helper objects
    my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;     
    my $ss = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;   
    my $is = Bio::Phylo::PhyLoTA::Service::InferenceService->new(
    	'tool'     => lc( $opt->inferencetool ),
    	'workdir'  => $opt->workdir,
    	'usertree' => $ts->make_usertree(
    		$supermatrix,
    		$starttree,
    		$self->workdir.'/user.dnd'
    	)
    );
    
    # run the analysis, process results
    my $base = $self->outfile;
    for my $i ( 1 .. $bootstrap ) {
    
    	# assign input matrix
        my $matrix;
        if ( $is->is_bayesian or $bootstrap == 1 ) {
        	$matrix = $supermatrix;
        }
        else {
         	$matrix = $ss->bootstrap( $supermatrix, $i );
        }
        
        # create replicate settings
		$self->outfile( "${base}.${i}" );
		$is->replicate( $i );
	    $is->configure;
	    
	    # run
		my $backbone = $is->run( 'matrix' => $matrix );  
		$self->_append( $backbone, $base . '.replicates' );
    }
    $self->_process_result( $base . '.replicates', !$opt->ids, $bootstrap - 1 );
    return 1;
}

sub _append {
    my ( $self, $in, $out ) = @_;
    $self->logger->info("appending contents of $in to $out");

    # read the in file
    open my $infh, '<', $in or die $!;
    my $content = do { local $/; <$infh> };
    close $infh;
    
    # write to the out file
    open my $outfh, '>>', $out or die $!;
    print $outfh $content;
    close $outfh;
}

# process the inferred tree, e.g. by mapping identifiers
sub _process_result {
    my ( $self, $backbone, $remap, $consense ) = @_;
    
    # build consensus if bootstrapping, otherwise read
    # single tree
    my $bbtree;
    if ( $consense ) {
        my $forest = parse(
            '-format' => 'newick',
            '-file'   => $backbone,
        );
        
        # here we will use treeannotator. XXX test to see if
        # TA will accept newick files.
        $bbtree = $forest->make_consensus( 
        	'-branches'  => 'average',
        	'-summarize' => 'fraction',
        );
    }
    else {
        $bbtree = parse_tree(
            '-format' => 'newick',
            '-file'   => $backbone,
        )
    }
    
    # map IDs to names, if requested
    if ( $remap ) {
        my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
        $bbtree = $ts->remap_to_name($bbtree); 
    }
    
    # write output file
    my $outfile = $self->outfile;
    open my $outfh, '>', $outfile or die $!;
    print $outfh $bbtree->to_newick( '-nodelabels' => 1 );
    close $outfh;   
    $self->logger->info("DONE, results written to $outfile");
}

1;
