#!/usr/bin/perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Bio::Phylo::IO qw(parse parse_tree);
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Tools::Run::Phylo::ExaML;
use Bio::Tools::Run::Phylo::ExaBayes;
use Bio::Phylo::PhyLoTA::Service::TreeService;

=head1 NAME

infer_backbone.pl - infers backbone tree using ExaML

=head1 SYNOPSIS

 $ infer_backbone.pl -o <file> -w <dir> -c <common tree> -s <supermatrix> -t <taxa> [--verbose]

=head1 DESCRIPTION

Given an input supermatrix (in interleaved PHYLIP format), infers a maximum likelihood tree 
using ExaML or ExaBayes and writes this to an output file. For tree inference with examl, an 
NCBI classification tree (in Newick format) has to be supplied as a starting tree.
ExaML and ExaBayes produce many intermediate checkpoint files, for which a
directory location needs to be specified.
Finally, the tree gets rerooted such that the number of monophyletic genera is maximized.

=cut

# process command line arguments
my $verbosity = WARN;
my ( $outfile, $workdir, $commontree, $supermatrix, $taxa );
GetOptions(
	'outfile=s'     => \$outfile,
	'workdir=s'     => \$workdir,
	'commontree=s'  => \$commontree,
	'supermatrix=s' => \$supermatrix,
        'taxa=s'        => \$taxa,
	'verbose+'      => \$verbosity,
);

# instantiate helper objects
my $config = Bio::Phylo::PhyLoTA::Config->new;
my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
my $logger = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => [qw(main Bio::Phylo::PhyLoTA::Service::TreeService 
                        Bio::Tools::Run::Phylo::ExaBayes)],
    );

# choose tool for backbone inference set in config file
my $tool;
if ( $config->BACKBONE_INFERENCE_TOOL eq "examl" ) {
        $logger->info("using examl for backbone inference");
        $tool  = Bio::Tools::Run::Phylo::ExaML->new;
} 
elsif ( $config->BACKBONE_INFERENCE_TOOL eq "exabayes" ) {
        $logger->info("using exabayes for backbone inference");
        $tool  = Bio::Tools::Run::Phylo::ExaBayes->new;
}
else {
        die ("tool for backbone inference must be set in config file.")
}

# set outfile name
$logger->info("going to create output file $outfile");
$tool->outfile_name($outfile);

# set working directory
$logger->info("going to use working directory $workdir");
$tool->work_dir($workdir);

# set parser location
$logger->info("going to use parser executable ".$config->PARSER_BIN);
$tool->parser($config->PARSER_BIN);

# set mpirun location
$logger->info("going to use mpirun executable ".$config->MPIRUN_BIN);
$tool->mpirun($config->MPIRUN_BIN);

# set number of nodes
$logger->info("setting number of MPI nodes ".$config->NODES);
$tool->nodes($config->NODES);

# set to quiet
$logger->info("setting backbone tool to requested verbosity");
$tool->quiet($verbosity <= WARN);


# set tool-specific arguments. if we infer with ExaML, the input tree is prepared
my $intree;

if ( $config->BACKBONE_INFERENCE_TOOL eq "examl" ) {
        # set examl location
        $logger->info("going to use executable ".$config->EXAML_BIN);
        $tool->executable($config->EXAML_BIN);
        
        # set substitution model
        $logger->info("setting substitution model ".$config->EXAML_MODEL);
        $tool->m($config->EXAML_MODEL);
   
        # here we need to read the names from the phylip file and adjust the
        # common tree accordingly: it must only retain the taxa in the supermatrix,
        # and it must be fully resolved.
        $intree = File::Spec->catfile( $workdir, 'user.dnd' );
        open my $fh, '>', $intree or die $!;
        print $fh parse(
                '-format' => 'newick',
                '-file'   => $commontree,
            )->first
            ->keep_tips( [ $ts->read_tipnames($supermatrix) ] )
            ->resolve
            ->remove_unbranched_internals
            ->deroot
            ->to_newick;
        close $fh;
}
elsif ( $config->BACKBONE_INFERENCE_TOOL eq "exabayes" ) {
        # set exabayes location
        $logger->info("going to use executable ".$config->EXABAYES_BIN);
        $tool->executable($config->EXABAYES_BIN);
        
        # set numbet of indepedent runs
        my $numruns = 4;
        $logger->info("going to perform $numruns independent runs");
        $tool->numRuns($numruns);

        # set number of coupled chains
        my $cupchains = 2;
        $logger->info("setting coupled chains to $cupchains");
        $tool->numCoupledChains($cupchains);
        
        # set number of generations
        my $numgen = 10;#00_000;
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
        $tool->consense_bin($config->CONSENSE_BIN);
        
        # set burnin fraction
        $logger->info("setting burnin fraction to ".$config->BURNIN);
        $tool->burnin_fraction($config->BURNIN);
}

$logger->info("patience please, running $supermatrix with starting tree");
if ( $intree ) {
        $logger->info("using tree $intree as starting point");
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

$logger->info("rerooting backbone tree such that the numebr of monophyletic genera is maximized");

my @records = $mt->parse_taxa_file($taxa);
my $rerooted = $ts->reroot_tree($bbtree, @records);

# write rerooted tree to specified output file
open my $outfh, '>', $outfile or die $!;
print $outfh $rerooted->to_newick;
close $outfh;

$logger->info("done, result written to $outfile");
