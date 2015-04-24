package Bio::SUPERSMART::App::smrt::Command::BBinfer::examl;
use strict;
use warnings;
use Bio::Tools::Run::Phylo::ExaML;
use Bio::SUPERSMART::App::smrt qw(-ignore);
use base 'Bio::SUPERSMART::App::smrt::Command::BBinfer';

=head1 NAME

examl.pm - wrapper for ExaML. No serviceable parts inside.

=cut

sub _configure {
    my ( $self, $tool, $config ) = @_;
    my $logger = $self->logger;
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

sub _create {
    my $self    = shift;
    my $logger  = $self->logger; 
    my $workdir = $self->workdir;
    my $tool    = Bio::Tools::Run::Phylo::ExaML->new;    

    # set working directory
    $logger->info("going to use working directory $workdir");
    $tool->work_dir($workdir);

    return $tool;
}

# infers a backbone tree with examl and returns the tree file. args:
# 'tool'   => wrapper,
# 'matrix' => alignment file name
# 'tree'   => user starting tree file name
sub _run {
    my ( $self, %args ) = @_;
    my $logger = $self->logger;
    my $backbone = $args{'tool'}->run(
            '-phylip' => $args{'matrix'},
            '-intree' => $self->_make_usertree(@args{qw[matrix tree]}),
    );
    $logger->info("examl tree inference was succesful") if $backbone;
    return $backbone;
}

1;