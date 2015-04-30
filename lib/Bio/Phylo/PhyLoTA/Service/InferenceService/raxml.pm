package Bio::Phylo::PhyLoTA::Service::InferenceService::raxml;
use strict;
use warnings;
use File::Spec;
use Bio::Tools::Run::Phylo::Raxml;
use Bio::Phylo::PhyLoTA::Service::InferenceService;
use base 'Bio::Phylo::PhyLoTA::Service::InferenceService';

=head1 NAME

Bio::Phylo::PhyLoTA::Service::InferenceService::raxml - Infers phylogenetic trees

=head1 DESCRIPTION

Provides functionality inferring trees and sets of trees by wrapping RAxML.

=over

=item create

Instantiates and configures the wrapper object, returns the instantiated wrapper.

=cut

sub create {
    my $self = shift;
    my $tool = Bio::Tools::Run::Phylo::Raxml->new;
    
    # configure raxml runner           
    $tool->w( $self->workdir );
    return $tool;
}

=item configure

Provides the $config object to the wrapper so that settings defined in $config object can 
be applied to the wrapper.

=cut

sub configure {
    my ( $self ) = @_;
    my $tool   = $self->wrapper;
    my $config = $self->config;
    $tool->outfile_name($self->outfile);      
    $tool->m($config->RAXML_MODEL);
    $tool->N($config->RAXML_RUNS);
    $tool->p($config->RANDOM_SEED);
    $tool->T($config->NODES);
}

=item run

Runs the analysis on the provided supermatrix. Returns a tree file.

=cut

sub run {
    my ( $self, %args ) = @_;
    my $t      = $self->wrapper;
    my $logger = $self->logger;
    
    # create output file name
    my $treefile = File::Spec->catfile( $t->w, 'RAxML_bestTree.' . $t->outfile_name );

    # configure bootstrap options   
    if ( $args{'boot'} ) {
        
        # set rapid bootstrap analysis
        $t->f('a');     
        
        # bootstrap random seed 
        $t->x($self->config->RANDOM_SEED);            
        $treefile = File::Spec->catfile( $t->w, 'RAxML_bipartitions.' . $t->outfile_name );     
    }   

    # run raxml, returns bioperl tree
    my $bptree = $t->run($args{'matrix'});          
    $logger->fatal('RAxML inference failed; outfile not present') if not -e $treefile; 
    return $treefile;
}

=back 

=cut

1;
