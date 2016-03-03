package Bio::SUPERSMART::App::smrtutils::Command::PruneTree;

use strict;
use warnings;

use Bio::SUPERSMART::Service::TreeService;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;

use Bio::Phylo::IO qw(parse_tree);

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: prunes user provided taxa from a tree

=head1 NAME

PruneTree.pm - prunes taxa from a given tree

=head1 SYNOPSYS

smrt-utils prunetree [-h ] [-v ] [-w <dir>] [-l <file>] [-y ] -t <file> [-o <file>] [-g ]

=head1 DESCRIPTION

Given an input tree and one or more taxa, prunes out the taxa from the tree. Input taxa of higher level are
expanded and all species (or lower) belonging to the higher level taxon are pruned out (if present in the tree).

=cut

sub options {

	my ($self, $opt, $args) = @_;
	my $outfile_default = 'pruned.dnd';
	my $format_default = 'nexus';
	return (
		['tree|t=s', 'file name of input tree', { arg => 'file', mandatory => 1}],
		['tree_format|f=s', "format of tree input file (newick or nexus), defaults to $format_default", { default => $format_default } ],
		['outfile|o=s', "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => 'file'}],
		['taxa|g=s', 'one or multiple taxa to be pruned from the tree', {} ],
		['negative_branches|n', 'prune out branches (and appendent nodes) which have negative length']
	    );
}

sub validate {
	my ($self, $opt, $args) = @_;

	# We only have to check the 'infile' argument.
	#  If the infile is absent or empty, abort
	my $file = $opt->tree;
	$self->usage_error('no tree argument given') if not $file;
	$self->usage_error('file $file does not exist') unless (-e $file);
	$self->usage_error('file $file is empty') unless (-s $file);
	$self->usage_error('taxa and/or negative_branches flag must be given') if not ($opt->taxa or $opt->negative_branches);
}

sub run {
	my ($self, $opt, $args) = @_;
	my $treefile = $opt->tree;
	my $taxa = $opt->taxa;
	my $outfile = $opt->outfile;
	my $logger = $self->logger;

	my $ts = Bio::SUPERSMART::Service::TreeService->new;
	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;

	my $tree = parse_tree(
		'-file'   => $treefile,
		'-format' => $opt->tree_format,
	    );
	my $ntax_initial = $tree->get_ntax;

	if ( $taxa ) {
		my @names = split(',', $taxa);
		for my $n (@names) {
			# taxa to prune could be higher level, so search all child taxa until lowest rank
			my %all_taxa = map {$_=~s/_/|/g; $_=~s/\ /_/g; $_=~s/\'|\"//g;  $_=>1} $mts->expand_taxa([$n], 'Forma');

			# select the taxa which are present in the tree
			my @to_prune;
			for my $node ( @{$tree->get_terminals} ) {
				my $terminal_name = $node->get_name;
				$terminal_name =~s/\'|\"//g;
				if ( $all_taxa{$terminal_name} ) {
					push @to_prune, $node->get_name;
				}
			}
			# prune
			$logger->info('Pruning ' . scalar (@to_prune) . " terminals from tree for taxon $n");
			$tree->prune_tips(\@to_prune);
		}
	}
	if ( $opt->negative_branches ) {
		my $count = 0;
		$tree->visit(sub{
			my $node = shift;
			my $bl = $node->get_branch_length || 0;
			if ( $bl < 0 ) {
				my @terminals = map {$_->get_name} @{$node->get_terminals};
				$count += scalar @terminals;
				$logger->info("Found negative branch length $bl, following tips will be removed: " . join(', ', @terminals));
				 $node->get_parent->prune_child( $node );
			}
			     });
		$tree->remove_unbranched_internals;
		$logger->info("Removed $count terminals while pruning negative branches");
	}

	my $ntax_pruned = $tree->get_ntax;
	$logger->info("Taxon number changed from $ntax_initial to $ntax_pruned");

	# write pruned tree to output file
	open my $out, '>', $outfile or die $!;
	print $out $tree->to_newick(nodelabels=>1);
	close $out;

	$logger->info("DONE. Outfile written to $outfile");
	return 1;
}

1;
