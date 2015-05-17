use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;

my $backbone = "${Bin}/testdata/piperaceae-backbone.dnd";
my $taxafile = "${Bin}/testdata/piperaceae-species.tsv";
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => ERROR,
	'-style' => 'simple',
	'-class' => [
		'main',
		'Bio::Phylo::PhyLoTA::Service::TreeService'
	]
);
	
# instantiate helper objects
my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;

# prepare taxa data
my @records = $mt->parse_taxa_file($taxafile);
my @ranks = ('forma', 'varietas', 'subspecies', 'species');	
my $level = $mt->get_highest_informative_level(@records);	
$log->info("highest informative taxon level : $level");

# iterate over trees
my $counter = 1;
open my $in, '<', $backbone or die $!;
while(<$in>) {
	my $newick_in = $_;
	my $tree = parse_tree( '-format' => 'newick', '-string' => $newick_in );	

	# use taxon IDs instead of names
	$ts->remap_to_ti($tree);
	
	# Perform rerooting on taxonomy			
	$tree->resolve;	
	$tree = $ts->reroot_tree($tree, \@records, [$level]);
	
	# smooth the basal branch
	$ts->smooth_basal_split($tree);
	
	# map back
	$tree = $ts->remap_to_name($tree);
	$tree->get_root->set_branch_length( 0.00 );
	
	# do the test
	my @negatives = grep { $_->get_branch_length < 0 } @{ $tree->get_entities };
	if ( @negatives ) {
		my $newick_out = $tree->to_newick( 'nodelabels' => 1 );
		ok( !@negatives, "have negative branches:\nin:\t$newick_in\nout:\t$newick_out" );
	}
	else {
		ok( 'no negative branches' );
	}
}