package Bio::SUPERSMART::App::smrtutils::Command::TreeDist;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse);
use Bio::Phylo::PhyLoTA::Service::TreeService; 

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: calculates the distance between two phylogenetic trees

=head1 NAME

Treedist.pm - Calculates the distance between two phylogenetic trees

=head1 SYNOPSYS

smrt-utils treedist [-t <file>] [-u <file>] [-f <format>] [-d] [-h] [-v] [-w <dir>] [-l <file>] [-y] 

=head1 DESCRIPTION

Calculates the distance between two phylogenetic trees. Currently supported distance measures: Robinson Foulds symmetric distance.

=cut

sub options {    
	my ($self, $opt, $args) = @_;
	my $format_default = 'newick';
	my $measure_default = 'rf';
	return (
		['tree1|t=s', "tree file", { arg => 'file' }],		
		['tree2|u=s', "tree file to compare to", { arg => 'file' }],		
		['treeformat|f=s', "file format of both input trees, defaults to $format_default", { default => $format_default, arg => "format" }],
		['distance_measure|d=s', "distance measure, currntly supported: rf (Robinson Foulds), defaults to $measure_default", { default => $measure_default}],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;			
	$self->usage_error('need two tree files as argument') if not ($opt->tree1 and $opt->tree2);
	$self->usage_error('tree file(s) not found') if not (-e $opt->tree1 and -e $opt->tree2);
	$self->usage_error('tree file(s) empty') if not (-s $opt->tree1 and -s $opt->tree2);
	
	my %measures = ( 'rf' => 1 );
	$self->usage_error('distance measure ' . $opt->distance_measure . ' not supported') if not ($measures{$opt->distance_measure});
}

sub run {
	my ($self, $opt, $args) = @_;    

	my $logger = $self->logger;      	
	
	my $tree1 = parse(
		'-file'   => $opt->tree1,
		'-format' => $opt->treeformat,
	    )->first;
	
	my $tree2 = parse(
		'-file'   => $opt->tree2,
		'-format' => $opt->treeformat,
	    )->first;

	# prune tips such that both trees have the same taxa.
	$logger->info("Pruning tips to have the same taxa in both trees");
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new; 
	my @trees = $ts->intersect_trees($tree1, $tree2);
	$tree1 = $trees[0];
	$tree2 = $trees[1];

	# calculate robinson-foulds distance and a normalized robinson foulds-distance
	my $diff;		
	my $ndiff;
	if ( $opt->distance_measure eq 'rf' ) { 
		$diff = $tree1->calc_symdiff($tree2);
		$ndiff = $tree1->calc_symdiff($tree2, 1);
	}
	
	$logger->info("Distance between " . $opt->tree1 . " and " . $opt->tree2 . " : $diff, normalized by splits: $ndiff");	
	$logger->info("DONE");
	
}



1;
