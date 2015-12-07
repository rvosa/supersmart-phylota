package Bio::SUPERSMART::Service::InferenceService::phyml;
use strict;
use warnings;
use Cwd 'abs_path';
use File::Spec;
use Bio::TreeIO;
use Bio::Tools::Run::Phylo::Phyml;
use Bio::SUPERSMART::Service::InferenceService;
use base 'Bio::SUPERSMART::Service::InferenceService';

=head1 NAME

Bio::SUPERSMART::Service::InferenceService::phyml - Infers phylogenetic trees

=head1 DESCRIPTION

Provides functionality inferring trees and sets of trees by wrapping PhyML.

=over

=item create

Instantiates and configures the wrapper object, returns the instantiated wrapper.

=cut

sub create {
    my $self = shift;
	return Bio::Tools::Run::Phylo::Phyml->new;
}

=item configure

Provides the $config object to the wrapper so that settings defined in $config object can 
be applied to the wrapper.

=cut

sub configure {
    my ( $self ) = @_;
    my $config = $self->config;
    my $tool   = $self->wrapper;
	$tool->data_type('nt');
	$tool->model($config->PHYML_MODEL);
	$tool->tree("BIONJ");
}

=item run

Runs the analysis on the provided supermatrix. Returns a tree file.

=cut

sub run {
    my ( $self, %args ) = @_;
    my $t      = $self->wrapper;
    my $logger = $self->logger;
	
    # run phyml, returns bioperl tree	
	$logger->info('Running phyml with matrix ' . $args{'matrix'});
    my $bptree = $t->run($args{'matrix'});
	
	$logger->fatal('No tree received from phyml') if not $bptree;

	# Write bioperl tree to file, 
	# file name is in the style of phyml output.
	my $treefile = $args{'matrix'} . "_phyml_tree.txt";
	my $out = new Bio::TreeIO(-file => ">$treefile",
							  -format => 'newick');	
	$out->write_tree($bptree);
	print {$out->_fh}  "\n";
    return $treefile;
}

=item cleanup

Cleans up any intermediate files.

=cut

sub cleanup { 
    my $self = shift;
    my $dir = $self->workdir;
    unlink $dir . '/*_phyml_tree.txt';
}

=back 

=cut

1;
