package Bio::Apps::Supersmart::Command::Classify;

use strict;
use warnings;

use Bio::Phylo::IO qw(unparse);
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;

use base 'Bio::Apps::GlobalCmd';
use Bio::Apps::Supersmart qw(-command);

# ABSTRACT: writes a classification tree for a given taxa table

=head1 NAME

Classify.pm - writes polytomous taxonomy tree

=head1 SYNOPSYS

=head1 DESCRIPTION

Given an input table of resolved taxa, produces a taxonomy tree, i.e. a polytomous 
cladogram with labeled interior nodes, in Newick format. This tree is used as a 
starting tree to infer a backbone phylogeny.

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "classification-tree.dnd";
	return (
		["infile|i=s", "file (tab-seperated value format) as produced by 'smrt taxize'", { required => 1, arg => "file"}],
		["outfile|o=s", "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],
		["nodelabels|n", "write labels for internal nodes", {default => 1}],
		["outformat|f=s", "format of generated tree, defaults to 'newick'", {default => 'newick'}],
	
	);	
}

sub validate_args {
	my ($self, $opt, $args) = @_;		

	# We only have to check the 'infile' argument. 
	#  If the infile is absent or empty, abort  
	my $file = $opt->infile;
	$self->usage_error("no infile argument given") if not $file;
	$self->usage_error("file $file does not exist") unless (-e $file);
	$self->usage_error("file $file is empty") unless (-s $file);			
}

sub execute {
	my ($self, $opt, $args) = @_;

	my $verbosity = INFO;
	
	# collect command-line arguments
	my $infile = $opt->infile;
	my $outfile = $opt->outfile;
	$verbosity += $opt->verbose ? $opt->verbose : 0;
	my $nodelabels = $opt->nodelabels;
	my $outformat = $opt->outformat;
	
	# instantiate helper objects
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $mt =  Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;

	my $log = Bio::Phylo::Util::Logger->new(
		'-class' => [ 'main', 'Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector',
					   'Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa' ],		
		'-level' => $verbosity
	);	
	# parse the taxa file 
	my @taxatable = $mt->parse_taxa_file($infile);

	# instantiate nodes from infile
	my @nodes = $mts->get_nodes_for_table( @taxatable );

	# compute common tree
	my $tree = $mts->get_tree_for_nodes(@nodes);
	$log->debug("done computing common tree");
	
	my $newick = $tree->to_newick;
	
	# create node labels
	$tree->visit(sub{
		my $node = shift;
		if ( $nodelabels or $node->is_terminal ) {
			my $label = $node->get_guid;		
			$node->set_name( $label );
		}
	});
	
	# write output
	open my $out, '>', $outfile or die $!;
	print $out unparse(
		'-format'     => $outformat,
		'-phylo'      => $tree,
		'-nodelabels' => $nodelabels,
	);
	close $out;
	$log->info("done, results written to $outfile");
	
	return 1;
}


1;