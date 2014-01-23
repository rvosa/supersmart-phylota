#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

# process command line arguments
my $verbosity = WARN;
my @levels = qw[species genus family order class phylum kingdom];
my ( $infile );
GetOptions(
	'infile=s' => \$infile,
	'verbose+' => \$verbosity,
	'level=s'  => \@levels,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => [qw(
		main
		Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector
	)]
);
my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;

# read names from file or STDIN, clean line breaks
my @names;
if ( $infile eq '-' ) {
	@names = <STDIN>;
	chomp(@names);
	$log->info("read species names from STDIN");
}
else {
	open my $fh, '<', $infile or die $!;
	@names = <$fh>;
	chomp(@names);
	$log->info("read ".scalar(@names)." species names from $infile");
}

# print table header
print join("\t", 'name', @levels), "\n";

# this will take some time to do the taxonomic name resolution in the
# database and with webservices
for my $name ( @names ) {
	my @nodes = $mts->get_nodes_for_names($name);
	if ( @nodes ) {
		if ( @nodes > 1 ) {
			$log->warn("Found more than one taxon for name $name");
		}
		
		# for each node, fetch the IDs of all taxonomic levels of interest
		for my $node ( @nodes ) {
			
			# create hash of taxonomic levels so that when we walk up the
			# taxonomy tree we can more easily check to see if we are at a
			# level of interest
			my %level = map { $_ => undef } @levels;
			
			# traverse up the tree
			while ( $node ) {
				my $rank = $node->rank;
				if ( exists $level{$rank} ) {
					$level{$rank} = $node->get_id;
				}
				$node = $node->get_parent;
			}
			print join("\t", $name, @level{@levels} ), "\n";
		}
	}
	else {
		$log->warn("Couldn't resolve name $name");
	}
}