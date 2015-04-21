package Bio::SUPERSMART::App::smrt::Command::BBinfer::exabayes;
use Bio::Tools::Run::Phylo::ExaBayes;
use base 'Bio::SUPERSMART::App::smrt::Command::BBinfer';

=head1 NAME

exabayes.pm - wrapper for ExaBayes. No serviceable parts inside.

=cut


sub _configure {
    my ( $self, $tool, $config ) = @_;
    my $logger = $self->logger;

    # set mpirun location
    $logger->info("going to use mpirun executable ".$config->MPIRUN_BIN);
    $tool->mpirun($config->MPIRUN_BIN);
    
    # set number of nodes
    $logger->info("setting number of MPI nodes ".$config->NODES);
    $tool->nodes($config->NODES);

    # set exabayes location
    $logger->info("going to use executable ".$config->EXABAYES_BIN);
    $tool->executable($config->EXABAYES_BIN);
        
    # set parser location
    $logger->info("going to use parser executable ".$config->EXABAYES_PARSER_BIN);
    $tool->parser($config->EXABAYES_PARSER_BIN);
        
    # set number of independent runs
    $logger->info("going to perform ".$config->EXABAYES_NUMRUNS." independent runs");
    $tool->numRuns($config->EXABAYES_NUMRUNS);

    # set number of coupled chains
    $logger->info("setting coupled chains to ".$config->EXABAYES_NUMCHAINS);
    $tool->numCoupledChains($config->EXABAYES_NUMCHAINS);
        
    # set number of generations
    $logger->info("setting number of generations to ".$config->EXABAYES_NUMGENS);
    $tool->numGen($config->EXABAYES_NUMGENS);
            
    # set random seed to config value
    $tool->seed($config->RANDOM_SEED); 
    
    # set binary of consense program
    $logger->info("setting consense binary");
    $tool->consense_bin($config->EXABAYES_CONSENSE_BIN);
        
    # set burnin fraction
    $logger->info("setting burnin fraction to ".$config->BURNIN);
    $tool->burnin_fraction($config->BURNIN); 
}

sub _create {
    my $self = shift;

    # instantiate helper objects and variables
    my $tool    = Bio::Tools::Run::Phylo::ExaBayes->new;    
    my $logger  = $self->logger;
    my $workdir = $self->workdir;
    my $outfile = $self->outfile;
        
    # set outfile name
    $logger->info("going to create output file $outfile");
    $tool->outfile_name($outfile);
    
    # set working directory
    $logger->info("going to use working directory $workdir");
    $tool->work_dir($workdir);        
        
    # set parsimony start
    $logger->info("setting ExaBayes to start from parsimony tree");
    $tool->parsimonyStart('true');
        
    # set run id to current process id
    $tool->run_id("infer_backbone_".$$);
        
    # set oufile format
    $logger->info("setting outfile format to 'newick'");
    $tool->outfile_format("newick");
        
    # set mode to zero (faster, but takes more memory)
    $logger->info("setting exabayes 'mode' argument to zero ");
    $tool->mode(0);             
    return $tool;
}

sub _run {
    my ( $self, %args ) = @_;
    return $args{'tool'}->run( '-phylip' => $args{'matrix'} );
}

sub _is_bayesian { 1 }

1;
