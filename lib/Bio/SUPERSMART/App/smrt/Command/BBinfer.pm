package Bio::SUPERSMART::App::smrt::Command::BBinfer;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse parse_tree);
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::TreeService;

use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: inference of a genus-level backbone tree

=head1 NAME

BBinfer.pm - inference of genus-level backbone tree

=head1 SYNOPSIS


=head1 DESCRIPTION

Given an input supermatrix (in interleaved PHYLIP format), infers either a maximum likelihood tree 
using ExaML or uses ExaBayes to sample a posterior distribution of trees. If ExaBayes is used, a 
'consense' tree is calculated from the posterior. The tree resulting from this script is written to file.
For tree inference with examl, an NCBI classification tree (in Newick format) has to be supplied 
as a starting tree. ExaML and ExaBayes produce many intermediate checkpoint files, for which a
directory location needs to be specified. The software for tree inference (ExaML or ExaBayes) is
determined in the configuration file, but can optionally also be given as a command-line argumet.

=cut


sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "backbone.dnd";
	return (
		["supermatrix|s=s", "matrix of concatenated multiple sequece alignments as produced by 'smrt bbmerge'", { arg => "file", mandatory => 1}],	
		["starttree|t=s", "starting tree for ML tree inference as for instance produced by 'smrt classify'", { arg => "file", mandatory => 1}],
		["inferencetool|i=s", "software tool for backbone inference, defaults to ExaML", {default => 'examl', arg => "tool"}],			
		["outfile|o=s", "name of the output tree file (in newick format), defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],			

	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	#  If alignment or taxa file is absent or empty, abort  
	my @files = ( $opt->supermatrix, $opt->starttree );
	foreach my $file ( @files ){
		$self->usage_error("need supermatrix and starttree arguments") if not $file;
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
}


sub run {
	my ($self, $opt, $args) = @_;		
		
	# collect command-line arguments
	my $supermatrix = $opt->supermatrix;
	my $outfile = $self->outfile;
	my $starttree = $opt->starttree;
	my $inferencetool = $opt->inferencetool;
	(my $workdir = $opt->workdir) =~ s/\/$//g;
	
	# instantiate helper objects
	my $config = Bio::Phylo::PhyLoTA::Config->new;
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $logger = $self->logger;
	
	# adjust capitalization if necessary
	$inferencetool = lc $inferencetool;
	$inferencetool =~ s/examl/ExaML/g;
	$inferencetool =~ s/exabayes/ExaBayes/g;
	
	my $toolclass = "Bio::Tools::Run::Phylo::" . $inferencetool;
	eval "use $toolclass";
	my $tool = $toolclass->new;
	
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
	
	# set tool-specific arguments. if we infer with ExaML, the input tree is prepared
	my $intree;

	if ( lc $inferencetool eq "examl" ) {
        # set examl location
        $logger->info("going to use executable ".$config->EXAML_BIN);
        $tool->executable($config->EXAML_BIN);

        # set parser location
        $logger->info("going to use parser executable ".$config->EXAML_PARSER_BIN);
        $tool->parser($config->EXAML_PARSER_BIN);

        # set substitution model
        $logger->info("setting substitution model ".$config->EXAML_MODEL);
        $tool->m($config->EXAML_MODEL);

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
		die("Number of taxa in supermatrix (" . scalar(@tipnames) .  
			") does not match number of taxa in NCBI classification tree ("
			 . scalar(@terminals) . ")") if not scalar(@terminals) == scalar(@tipnames);

		$intree = File::Spec->catfile( $workdir, 'user.dnd' );
        open my $fh, '>', $intree or die $!;
        print $fh $tree->to_newick;
		close $fh;
	
	}
elsif ( lc $inferencetool eq "exabayes" ) {
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
        my $numgen = 1_000_000;
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
	}
	$logger->info("patience please, running backbone inference with $supermatrix");
	
	if ( $intree ) {
	        $logger->info("using tree $intree as starting tree");
	}
	
	my $backbone = $tool->run(
	        '-phylip' => $supermatrix,
	        '-intree' => $intree,
	    );
	
	
	# read generated backbone tree from file ane reroot it
	my $bbtree = parse_tree(
		'-format'     => 'newick',
		'-file'       => $backbone,
		'-as_project' => 1,
	);
	
	$ts->remap_to_name($bbtree);
	
	open my $outfh, '>', $outfile or die $!;
	print $outfh $bbtree->to_newick;
	close $outfh;
	
	
	
	$logger->info("DONE, results written to $outfile");
	
}

1;