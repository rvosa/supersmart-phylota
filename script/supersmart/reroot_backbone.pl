#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Bio::Phylo::IO qw(parse parse_tree);

use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::TreeService;

=head1 NAME

reroot_backbone.pl - reroots backbone tree such that the amount of paraphyletic species is minimized

=head1 SYNOPSIS

 $ reroot_backbone.pl -b <treefile> -w <dir> -t <taxa> [--verbose]

=head1 DESCRIPTION

Given an input backbone tree in newick format and a list of taxa with their NCBI taxonomy identifiers
for different taxonomic ranks, re-roots the tree such that the amount of paraphyletic species
(with respect to certain taxonomic ranks) is minimized. Output is the rerooted tree in newick format.

=cut


# process command line arguments
my $verbosity = INFO;
my ( $backbone, $workdir, $taxa );
GetOptions(
	'backbone=s'    => \$backbone,
	'workdir=s'     => \$workdir,
    'taxa=s'        => \$taxa,
	'verbose+'      => \$verbosity,
);

# instantiate helper objects
my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
my $logger = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => [ qw(main Bio::Phylo::PhyLoTA::Service::TreeService) ], 
	);

my $tree = parse_tree(
	'-file'   => $backbone,
	'-format' => 'newick',
);

my @records = $mt->parse_taxa_file($taxa);

# reroot the tree
$logger->info("rerooting backone tree");
my $rerooted = $ts->reroot_tree($tree, @records);

#print the result
print $rerooted->to_newick;

$logger->info("done rerooting tree");
