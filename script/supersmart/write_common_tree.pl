#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'unparse';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

# process command line arguments
my $template   = '${get_guid}';  # by default, the node label is the NCBI taxon ID
my @properties = qw(get_guid); # additional column values to use in template
my $infile     = '-'; # read from STDIN by default
my $outformat  = 'newick'; # write newick by default
my $nodelabels = 1; # write internal node labels by default
my $verbosity  = WARN; # low verbosity
GetOptions(
	'template=s'  => \$template,
	'infile=s'    => \$infile,
	'outformat=s' => \$outformat,
	'verbose+'    => \$verbosity,
	'property=s'  => \@properties,
	'nodelabels'  => $nodelabels,
);

# instantiate helper objects
my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-class' => [ qw(main Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector) ],
	'-level' => $verbosity
);

# instantiate nodes from infile
my @nodes = $mts->get_nodes_for_table( '-file' => $infile );

# compute common tree
my $tree = $mts->get_tree_for_nodes(@nodes);
$log->debug("done computing common tree");

# create node labels
$tree->visit(sub{
	my $node = shift;
	if ( $nodelabels or $node->is_terminal ) {
		my $label;		
		my $statement = "my (" . join( ',', map ("\$$_", @properties) ). ");\n";
		for my $property ( @properties ) {
			$statement .= "\$$property = q[" . $node->$property . "];\n";
			$log->debug($statement);
		}		
		$statement .= "\$label = $template";
		$log->debug($statement);
		eval $statement;
		die $@ if $@;
		$log->debug($label);
		$node->set_name( $label );
	}
});

# write output
print unparse(
	'-format'     => $outformat,
	'-phylo'      => $tree,
	'-nodelabels' => $nodelabels,
);