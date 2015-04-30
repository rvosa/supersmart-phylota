package Bio::Phylo::PhyLoTA::Service::InferenceService::exabayes;
use Bio::Tools::Run::Phylo::ExaBayes;
use Bio::Phylo::PhyLoTA::Service::InferenceService;
use base 'Bio::Phylo::PhyLoTA::Service::InferenceService';

=head1 NAME

Bio::Phylo::PhyLoTA::Service::InferenceService::exabayes - Infers phylogenetic trees

=head1 DESCRIPTION

Provides functionality inferring trees and sets of trees by wrapping ExaBayes.

=over

=item configure

Provides the $config object to the wrapper so that settings defined in $config object can 
be applied to the wrapper.

=cut

sub configure {
    my ( $self ) = @_;
    my $tool    = $self->wrapper;
    my $config  = $self->config;
    my $logger  = $self->logger;
    my $outfile = $self->outfile;
        
    # set outfile name
    $logger->info("going to create output file $outfile");
    $tool->outfile_name($outfile);    

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

=item create

Instantiates and configures the wrapper object, returns the instantiated wrapper.

=cut

sub create {
    my $self = shift;

    # instantiate helper objects and variables
    my $logger  = $self->logger;
    my $workdir = $self->workdir;
    my $tool    = Bio::Tools::Run::Phylo::ExaBayes->new;
    
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

=item run

Runs the analysis on the provided supermatrix. Returns a tree file.

=cut

sub run {
    my ( $self, %args ) = @_;
    return $self->wrapper->run( '-phylip' => $args{'matrix'} );
}

=item is_bayesian

Indicates whether the inference service is bayesian, in which case we will
not run a bootstrapping analysis. Returns true.

=cut

sub is_bayesian { 1 }

=item cleanup

Cleans up any intermediate files.

=cut

sub cleanup {
    my $self = shift;
    my $runid = $self->wrapper->run_id;
    my $dir = $self->workdir;    

    # remove the input matrices
    unlink "${dir}/${runid}-dat.binary";
    unlink "${dir}/${runid}.nex";

    # remove single intermediate files
    for my $prefix ( qw(diagnostics checkpoint info prevCheckpointBackup ConsensusExtendedMajorityRuleNexus) ) {
        unlink "${dir}/ExaBayes_${prefix}.${runid}";
    }

    # remove files from chains
    opendir my $dh, $dir or die $!;
    while( my $entry = readdir $dh ) {
        if ( $entry =~ /ExaBayes_parameters.${runid}.\d+/ ) {
            unlink "${dir}/${entry}";
        }
    }
}

=back 

=cut

1;
