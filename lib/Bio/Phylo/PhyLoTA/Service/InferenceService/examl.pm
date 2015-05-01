package Bio::Phylo::PhyLoTA::Service::InferenceService::examl;
use strict;
use warnings;
use Bio::Tools::Run::Phylo::ExaML;
use Bio::Phylo::PhyLoTA::Service::InferenceService;
use base 'Bio::Phylo::PhyLoTA::Service::InferenceService';

=head1 NAME

Bio::Phylo::PhyLoTA::Service::InferenceService::examl - Infers phylogenetic trees

=head1 DESCRIPTION

Provides functionality inferring trees and sets of trees by wrapping ExaML.

=over

=item configure

Provides the $config object to the wrapper so that settings defined in $config object can 
be applied to the wrapper.

=cut

sub configure {
    my ( $self ) = @_;
    my $logger  = $self->logger;
    my $outfile = $self->outfile;   
    my $tool    = $self->wrapper;
    my $config  = $self->config; 
    
    # set outfile name
    $logger->info("going to create output file $outfile");
    $tool->outfile_name($outfile);    

    # set mpirun location
    $logger->info("going to use mpirun executable ".$config->MPIRUN_BIN);
    $tool->mpirun($config->MPIRUN_BIN);
    
    # set number of nodes
    $logger->info("setting number of MPI nodes ".$config->NODES);
    $tool->nodes($config->NODES);
    
    # set executable location
    $logger->info("going to use executable ".$config->EXAML_BIN);
    $tool->executable($config->EXAML_BIN);

    # set parser location
    $logger->info("going to use parser executable ".$config->EXAML_PARSER_BIN);
    $tool->parser($config->EXAML_PARSER_BIN);

    # set substitution model
    $logger->info("setting substitution model ".$config->EXAML_MODEL);
    $tool->m($config->EXAML_MODEL);

    # set random seed
    $tool->p($config->RANDOM_SEED);
}

=item create

Instantiates and configures the wrapper object, returns the instantiated wrapper.

=cut

sub create {
    my $self    = shift;
    my $logger  = $self->logger; 
    my $workdir = $self->workdir;
    my $tool    = Bio::Tools::Run::Phylo::ExaML->new;    

    # set working directory
    $logger->info("going to use working directory $workdir");
    $tool->work_dir($workdir);

    return $tool;
}

=item run

Runs the analysis on the provided supermatrix. Returns a tree file.

=cut

sub run {
    my ( $self, %args ) = @_;
    my $logger = $self->logger;
    $self->wrapper->run_id( 'examl-run-' . $$ . '-' . $self->replicate );
    my $backbone = $self->wrapper->run(
            '-phylip' => $args{'matrix'},
            '-intree' => $self->usertree,
    );
    $logger->info("examl tree inference was succesful") if $backbone;
    return $backbone;
}

=item cleanup

Cleans up all intermediate files.

=cut

sub cleanup { 
    my $self = shift;
    my ( $v, $d, $f ) = File::Spec->splitpath( $self->outfile );
    my $tempfile = 'ExaML_modelFile.' . $f;
    if ( -e $tempfile ) {
        unlink $tempfile;
        $self->logger->info("removed $tempfile");
    }
}

=back 

=cut

1;
