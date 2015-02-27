package Bio::SUPERSMART::App::smrt::Command::BBinfer;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse parse_tree);
use Bio::Phylo::Factory;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Tools::Run::Phylo::Raxml;
use Bio::Tools::Run::Phylo::ExaML;
use Bio::Tools::Run::Phylo::ExaBayes;
	
use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: inference of a genus-level backbone tree

=head1 NAME

BBinfer.pm - inference of genus-level backbone tree

=head1 SYNOPSIS

smrt bbinfer [-h ] [-v ] [-w <dir>] -s <file> -t <file> [-i <tool>] [-o <file>] 

=head1 DESCRIPTION

Given an input supermatrix (in interleaved PHYLIP format), infers either a maximum likelihood tree 
using ExaML or RaXML or uses ExaBayes to sample a posterior distribution of trees. If ExaBayes is used, a 
'consense' tree is calculated from the posterior. The tree resulting from this script is written to file.
For tree inference with examl, an NCBI classification tree (in Newick format) has to be supplied 
as a starting tree. ExaML and ExaBayes produce many intermediate checkpoint files, for which a
directory location needs to be specified. The software for tree inference (ExaML or ExaBayes) is
determined in the configuration file, but can optionally also be given as a command-line argumet.

=cut


sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "backbone.dnd";
	my $tool_default = "examl";
	return (
		["supermatrix|s=s", "matrix of concatenated multiple sequece alignments as produced by 'smrt bbmerge'", { arg => "file", mandatory => 1}],	
		["starttree|t=s", "starting tree for ML tree inference as for instance produced by 'smrt classify', only needed when inferencetool is ExaML", { arg => "file", mandatory => 1}],
		["inferencetool|i=s", "software tool for backbone inference (RaXML, ExaML or ExaBayes), defaults to $tool_default", {default => $tool_default, arg => "tool"}],			
		["bootstrap|b", "do a bootstrap analysis and add the support values to the backbone tree. Currently only supported for inferencetool RaXML", {}],
		["outfile|o=s", "name of the output tree file (in newick format), defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],			

	);	
}

sub validate {
	my ($self, $opt, $args) = @_;

	my $sm = $opt->supermatrix;
	my $st = $opt->starttree;
	my $tool = $opt->inferencetool;
	
	$self->usage_error("need supermatrix") if not $sm;
	$self->usage_error("file $sm does not exist") unless (-e $sm);
	$self->usage_error("file $sm is empty") unless (-s $sm);
	
	if ( lc $tool eq 'examl') {
		$self->usage_error("need starting tree") if not $st;
		$self->usage_error("file $st does not exist") unless (-e $st);
		$self->usage_error("file $st is empty") unless (-s $st);		
	}
}


sub run {
	my ($self, $opt, $args) = @_;		
		
	# collect command-line arguments
	my $supermatrix = $opt->supermatrix;
	my $outfile = $self->outfile;
	my $starttree = $opt->starttree;
	my $inferencetool = $opt->inferencetool;
	
	# instantiate helper objects
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $logger = $self->logger;
	
	my $backbone;
	
	if ( lc $inferencetool eq "examl") {
		$backbone = $self->_infer_examl ($supermatrix, $starttree);
	}	
	elsif ( lc $inferencetool eq "raxml") {
		$backbone = $self->_infer_raxml ($supermatrix, $opt->bootstrap);
	}	
	elsif ( lc $inferencetool eq "exabayes") {
		$backbone = $self->_infer_exabayes ($supermatrix);
	}
	
	# read generated backbone tree from file
	my $bbtree = parse_tree(
		'-format'     => 'newick',
		'-file'       => $backbone,
		'-as_project' => 1,
	);
	
	$ts->remap_to_name($bbtree);
	
	open my $outfh, '>', $outfile or die $!;
	print $outfh $bbtree->to_newick('-nodelabels' => 1);
	#print $outfh $ts->write_newick_tree($bbtree);

	close $outfh;
	
	$logger->info("DONE, results written to $outfile");
	
	return 1;
}


# infers a backbone tree with raxml and returns the tree file
sub _infer_raxml {
	my ($self, $supermatrix, $bootstrap) = @_;
	
	my $config = Bio::Phylo::PhyLoTA::Config->new;
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	
	
	
	my $tool = Bio::Tools::Run::Phylo::Raxml->new(-N => 100, -p => 1, -T => $config->NODES	);			
	$tool->outfile_name("out");   			
	$tool->w( $tool->tempdir );
	$tool->m("GTRGAMMA");	
	
	if ( $bootstrap ) {
		# set rapid bootstrap analysis
		$tool->f('a');		
		# bootstrap random seed 
		$tool->x(1);			
	}
		
	my $bptree = $tool->run($supermatrix);		
	##$tool->cleanup;
	
	# bioperl gives back a tree object, 
	# and to keep consistent with the exabayes and examl subroutines,
	# we write the tree to the specified outfile
	my $tree = Bio::Phylo::Forest::Tree->new_from_bioperl( $bptree );
		
	open my $fh, '>', $self->outfile;
	print $fh $tree->to_newick('-nodelabels' => 1);	
	close $fh;
	
	return $self->outfile;	

}

# infers a backbone tree with examl and returns the tree file
sub _infer_examl {
	my ($self, $supermatrix, $starttree) = @_;

	my $config = Bio::Phylo::PhyLoTA::Config->new;

	my $workdir = $self->workdir;
	my $tool = Bio::Tools::Run::Phylo::ExaML->new;
	
	my $logger = $self->logger;
	my $outfile = $self->outfile;
	
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	
	# set outfile name
	$logger->info("going to create output file $outfile");
	$tool->outfile_name($outfile);
	
	# set working directory
	$logger->info("going to use working directory $workdir");
	$tool->work_dir($workdir);
	
	# set mpirun location
	$logger->info("going to use mpirun executable ".$config->MPIRUN_BIN);
	$tool->mpirun($config->MPIRUN_BIN);
	
	# set number of nodes
	$logger->info("setting number of MPI nodes ".$config->NODES);
	$tool->nodes($config->NODES);
	
	# set to quiet
	#$logger->info("setting backbone tool to requested verbosity");
	#$tool->quiet('-v');
	
	$logger->info("going to use executable ".$config->EXAML_BIN);
    $tool->executable($config->EXAML_BIN);

    # set parser location
	$logger->info("going to use parser executable ".$config->EXAML_PARSER_BIN);
    $tool->parser($config->EXAML_PARSER_BIN);

    # set substitution model
    $logger->info("setting substitution model ".$config->EXAML_MODEL);
    $tool->m($config->EXAML_MODEL);

	# set random seed
	$tool->p(1); 

    # here we need to read the names from the phylip file and adjust the
    # common tree accordingly: it must only retain the taxa in the supermatrix,
    # and it must be fully resolved.
		
	my @tipnames = $ts->read_tipnames($supermatrix);

	my $tree = parse(
                '-format' => 'newick',
                '-file'   => $starttree,
        )->first;

    $ts->remap_to_ti($tree);
		           
    $tree->keep_tips( \@tipnames )
            ->resolve
            ->remove_unbranched_internals
            ->deroot;		

	my @terminals = @{$tree->get_terminals()};
	
	# it can occur that a taxon which has been chosen as an exemplar is
	#  not in the classification tree. For example, if a species has one subspecies,
	#  but the species and not the subspecies is an exemplar. Then this node
	#  is an unbranched internal and not in the classification tree. We therefore
	#  add a node to the starting tree

	if ( not scalar(@terminals) == scalar(@tipnames) ) {
		$logger->warn("Number of taxa in supermatrix (" . scalar(@tipnames) .  
			") does not match number of taxa in starting tree ("
			 . scalar(@terminals) . "). Could one of the exemplars be an internal node?");

		# get taxa that are in the supermatrix but not in the classification tree
		
		my @terminal_ids = map{ $_->get_name } @terminals;	
		my %diff;
		@diff {@terminal_ids} = ();
		my @diffs =  grep !exists($diff{$_}), @tipnames;
		
		# insert new node into starting tree
		
		my $fac = Bio::Phylo::Factory->new;
		foreach my $d (@diffs){
			$logger->warn("Adding node $d (" . $ts->find_node($d)->get_name . ") to starting tree");
			my $node = $fac->create_node('-guid' => $d, '-name' => $d);
			$node->set_parent(@{$tree->get_internals}[0]);
			$tree->insert($node);
			$tree->resolve;
		}	
	}

	my $intree = File::Spec->catfile( $workdir, 'user.dnd' );
    open my $fh, '>', $intree or die $!;
    print $fh $tree->to_newick;
	close $fh;
	$logger->info("Writing starting tree for examl inference to $intree");
	
	my $backbone = $tool->run(
	        '-phylip' => $supermatrix,
	        '-intree' => $intree,
	 );
	 $logger->info("examl tree inference was succesful");
	return $backbone;
}

# infers a backbone tree with exabayes and returns the tree file
sub _infer_exabayes {
	my ($self, $supermatrix) = @_;

	my $config = Bio::Phylo::PhyLoTA::Config->new;

	my $workdir = $self->workdir;
	my $tool = Bio::Tools::Run::Phylo::ExaBayes->new;
	
	my $logger = $self->logger;
	my $outfile = $self->outfile;
	
	# set outfile name
	$logger->info("going to create output file $outfile");
	$tool->outfile_name($outfile);
	
	# set working directory
	$logger->info("going to use working directory $workdir");
	$tool->work_dir($workdir);
	
	# set mpirun location
	$logger->info("going to use mpirun executable ".$config->MPIRUN_BIN);
	$tool->mpirun($config->MPIRUN_BIN);
	
	# set number of nodes
	$logger->info("setting number of MPI nodes ".$config->NODES);
	$tool->nodes($config->NODES);
	
	# set to quiet
	#$logger->info("setting backbone tool to requested verbosity");
	#$tool->quiet('-v');

    # set exabayes location
    $logger->info("going to use executable ".$config->EXABAYES_BIN);
    $tool->executable($config->EXABAYES_BIN);
        
    # set parser location
    $logger->info("going to use parser executable ".$config->EXABAYES_PARSER_BIN);
    $tool->parser($config->EXABAYES_PARSER_BIN);
        
    # set numbet of indepedent runs
    my $numruns = 2;
    $logger->info("going to perform $numruns independent runs");
    $tool->numRuns($numruns);

    # set number of coupled chains
    my $cupchains = 2;
    $logger->info("setting coupled chains to $cupchains");
    $tool->numCoupledChains($cupchains);
        
    # set number of generations
    my $numgen = 1_00_000;
    $logger->info("setting number of generations to $numgen");
    $tool->numGen($numgen);
        
    # set parsimony start
    $logger->info("setting ExaBayes to start from parsimony tree");
    $tool->parsimonyStart('true');
        
    # set run id to current process id
    $tool->run_id("infer_backbone_".$$);
        
    # set random seed to 1
    $tool->seed(1);
        
    # set oufile format
    $logger->info("setting outfile format to 'newick'");
    $tool->outfile_format("newick");
        
    # set binary of consense program
    $logger->info("setting consense binary");
    $tool->consense_bin($config->EXABAYES_CONSENSE_BIN);
        
    # set burnin fraction
    $logger->info("setting burnin fraction to ".$config->BURNIN);
    $tool->burnin_fraction($config->BURNIN);

    # set mode to zero (faster, but takes more memory)
    $logger->info("setting exabayes 'mode' argument to zero ");
    $tool->mode(0);
        
    # try to save memory with SEV technique
    #$logger->info("trying to save memory with SEV technique in exabayes run");
        #$tool->S(1);                	

	my $backbone = $tool->run(
	        '-phylip' => $supermatrix,
	    );
		
	return $backbone;
}

1;