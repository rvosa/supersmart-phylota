#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

# process command line arguments
my ( $format, $verbosity, $infile ) = ( 'newick', WARN );
GetOptions(
	'infile=s' => \$infile,
	'format=s' => \$format,
	'verbose+' => \$verbosity,
);

# create helper object
my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# parse tree
my $tree = parse_tree(
	'-format' => $format,
	'-file'   => $infile,
	'-as_project' => 1,
);

# do the renaming
$tree->visit(sub{
	my $node = shift;
	if ( $node->is_terminal ) {
		if ( my $name = $node->get_name ) {
			if ( $name =~ /^t?(\d+)$/ ) {
				my $tid = $1;
				my $taxon = $mts->find_node($tid);
				my $binomial = $taxon->taxon_name;
				$log->info("found name $binomial for ID $tid");
				$binomial =~ s/ /_/g;				
				$node->set_name($binomial);
			}
			else {
				$log->warn("$name doesn't match expectation");
			}
		}
		else {
			$log->warn("$node has no name");
		}
	}
});

# print output
print $tree->to_newick;

