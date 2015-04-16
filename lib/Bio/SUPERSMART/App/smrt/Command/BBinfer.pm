package Bio::SUPERSMART::App::smrt::Command::BBinfer;
use strict;
use warnings;
use Bio::Phylo::IO qw(parse parse_tree);
use Bio::Phylo::Factory;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::Util::Exceptions 'throw';
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: inference of a genus-level backbone tree

=head1 NAME

BBinfer.pm - inference of genus-level backbone tree

=head1 SYNOPSIS

smrt bbinfer [-h ] [-v ] [-w <dir>] -s <file> -t <file> [-i <tool>] [-o <file>] 

=head1 DESCRIPTION

Given an input supermatrix (in interleaved PHYLIP format), infers either a maximum 
likelihood tree using ExaML or RAxML or uses ExaBayes to sample a posterior distribution 
of trees. If ExaBayes is used, a 'consense' tree is calculated from the posterior. The 
tree resulting from this script is written to file. For tree inference with examl, an NCBI 
classification tree (in Newick format) has to be supplied as a starting tree. ExaML and 
ExaBayes produce many intermediate checkpoint files, for which a directory location needs 
to be specified. 

=cut

# define additional command line arguments
sub options {
    my ($self, $opt, $args) = @_;       
    my $outfile_default = "backbone.dnd";
    my $tool_default = "examl";
    return (
        ["supermatrix|s=s", "matrix of concatenated multiple sequece alignments as produced by 'smrt bbmerge'", { arg => "file", mandatory => 1}],  
        ["starttree|t=s", "starting tree for ExaML tree inference. If not given, a random starting tree is generated", { arg => "file", mandatory => 1}],
        ["inferencetool|i=s", "software tool for backbone inference (RAxML, ExaML or ExaBayes), defaults to $tool_default", {default => $tool_default, arg => "tool"}],         
        ["bootstrap|b", "do a bootstrap analysis and add the support values to the backbone tree. Currently only supported for inferencetool RAxML", {}],
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
    $self->usage_error("file $sm does not exist") unless (-e $sm);
    $self->usage_error("file $sm is empty") unless (-s $sm);
}

# run the analysis
sub run {
    my ( $self, $opt, $args ) = @_;     
        
    # collect command-line arguments
    my $supermatrix = $opt->supermatrix;
    my $starttree   = $opt->starttree;
    my $bootstrap   = $opt->bootstrap;
    my $toolname    = lc $opt->inferencetool;
    
    # check to see if we should become an instance of a subclass instead
    if ( ref($self) !~ /$toolname$/ ) {
        my $subclass = __PACKAGE__ . '::' . $toolname;
        eval "require $subclass";
        $@ and throw 'ExtensionError' => "Can't load $toolname wrapper: ".$@;
        bless $self, $subclass;
    }
    
    # instantiate and configure wrapper object
    my $tool = $self->_create;
    $self->_configure( $tool, $self->config );
    
    # run the analysis, process results
    my $backbone = $self->_run(
        'tool'   => $tool,
        'matrix' => $supermatrix,
        'tree'   => $starttree,
        'boot'   => $bootstrap,
    );  
    $self->_process_result($backbone, !$opt->ids);
    return 1;
}

# instantiates and configures the wrapper object,
# returns the instantiated wrapper
sub _create {
    throw 'NotImplemented' => "missing _create method in child class " . ref(shift);
}

# provides the $config object to the wrapper so that settings defined in $config
# object can be applied to the wrapper
sub _configure {
    throw 'NotImplemented' => "missing _confgure method in child class " . ref(shift);
}

# runs the analysis on the provided supermatrix. can potentially be run multiple
# times (in parallel?) to do bootstrapping. returns a tree file.
sub _run {
    throw 'NotImplemented' => "missing _run method in child class " . ref(shift);
}

# process the inferred tree, e.g. by mapping identifiers
sub _process_result {
    my ( $self, $backbone, $remap ) = @_;
    
    # read generated backbone tree from file
    my $bbtree = parse_tree(
        '-format' => 'newick',
        '-file'   => $backbone,
    );
    
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

# reconcile the common tree with the taxa in the supermatrix and prepare it
# for usage (remove unbranched internal nodes, randomly resolve polytomies, deroot)
sub _process_commontree {
    my ( $self, $commontree, @tipnames ) = @_;
    my $logger = $self->logger;
    my $ts     = Bio::Phylo::PhyLoTA::Service::TreeService->new;

    # read the common tree, map names to IDs, 
    # only retain tips in supermatrix
    my $tree = parse_tree(
        '-format' => 'newick',
        '-file'   => $commontree,
    );
    $ts->remap_to_ti( $tree );
    $tree->keep_tips( \@tipnames );
            
    # it can occur that a taxon that has been chosen as an exemplar is
    # not in the classification tree. For example, if a species has one subspecies,
    # but the species and not the subspecies is an exemplar. Then this node
    # is an unbranched internal and not in the classification tree. We therefore
    # add these node(s) to the starting tree
    my @terminals = @{ $tree->get_terminals };
    if ( @terminals != @tipnames ) {
        $logger->warn("tips mismatch: ".scalar(@tipnames)."!=".scalar(@terminals));

        # insert unseen nodes into starting tree        
        my %seen = map { $_->get_name => 1 } @terminals;        
        my $fac  = Bio::Phylo::Factory->new;
        my ($p)  = @{ $tree->get_internals };
        for my $d ( grep { ! $seen{$_} } @tipnames ) {
            $logger->warn("Adding node $d (" . $ts->find_node($d)->get_name . ") to starting tree");
            my $node = $fac->create_node( '-name' => $d );
            $node->set_parent($p);
            $tree->insert($node);
        }
    }
        
    # finalize the tree
    $tree->resolve->remove_unbranched_internals->deroot;
    return $tree;
}

# make a starting tree either by simulation or from input tree
sub _make_usertree {
    my ($self, $supermatrix, $commontree) = @_;
    my $logger   = $self->logger;
    my $ts       = Bio::Phylo::PhyLoTA::Service::TreeService->new;
    my @tipnames = $ts->read_tipnames($supermatrix);
    
    # this will be the tree object to write to file         
    my $tree;
    
    # process the common tree to the right set of tips, resolve it,
    # remove unbranched internals (due to taxonomy) and remove the root
    if ( $commontree ) {
        $logger->info("Going to prepare starttree $commontree for usage");
        $tree = $self->_process_commontree($commontree,@tipnames);
    }
    
    # simulate a random tree
    else {
        $logger->info("No starttree given as argument. Generating random starting tree.");  
        $tree = $self->_make_random_starttree(\@tipnames);  
    }
    
    # write to file
    my $intree = File::Spec->catfile( $self->workdir, 'user.dnd' );
    open my $fh, '>', $intree or die $!;
    print $fh $tree->to_newick();
    return $intree; 
}

# simulate an equiprobable starting tree
sub _make_random_starttree {
    my ( $self, $tipnames) = @_;
    my $fac  = Bio::Phylo::Factory->new;
    my $taxa = $fac->create_taxa;
    my $tree = $fac->create_tree;
    my $rootnode = $fac->create_node( '-name' => 'root' );
    $tree->insert($rootnode);
    for my $t (@{$tipnames}) {
        my $node = $fac->create_node( '-name' => $t, '-branch_length' => 0 );
        $node->set_parent($rootnode);
        $tree->insert($node);
    }
    $tree->resolve;
    $self->logger->info($tree->to_newick);
    return $tree;
}

1;
