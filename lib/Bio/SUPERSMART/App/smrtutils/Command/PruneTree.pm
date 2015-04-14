package Bio::SUPERSMART::App::smrtutils::Command::PruneTree;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

use Bio::Phylo::IO qw(parse);

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
	return (
		['tree|t=s', 'file name of input tree (newick format)', { arg => 'file', mandatory => 1}],
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
	
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	
	my $tree = parse(
		'-file'   => $treefile,
		'-format' => 'newick',
	    )->first;
	
	my $ntax_initial = $tree->get_ntax;
	
	if ( $taxa ) {		    
		my @names = split(',', $taxa);    
		for my $n (@names) {
			# taxa to prune could be higher level, so search all child taxa until lowest rank
			my %all_taxa = map {$_=~s/_/|/g; $_=~s/\ /_/g; $_=>1} $mts->expand_taxa([$n], 'Forma');
			#print "Found " . scalar(keys(%all_taxa)) . " taxa \n";
			#use Data::Dumper;
			#print Dumper(\%all_taxa);
			# select the taxa which are present in the tree
			my @nodes = grep {exists ($all_taxa{$_->get_name}) } @{$tree->get_terminals};				    
			#print "Number of nodes : " . scalar(@nodes) . "\n";
			my @to_prune = map {$_->get_name} @nodes; 				    
			
			# prune
			$logger->info('Pruning ' . scalar (@to_prune) . " terminals from tree for taxon $n");				    
			$tree->prune_tips(\@to_prune);
		}
	}      
	if ( $opt->negative_branches ) {
		$tree->visit(sub{ 
			my $node = shift;
			my $bl = $node->get_branch_length || 0;
			if ( $bl < 0 ) {
				my @terminals = map {$_->get_name} @{$node->get_terminals};
				$logger->info("Found negative branch length $bl, following tips will be removed: " . join(', ', @terminals));
				 $node->get_parent->prune_child( $node );
			}
			     });
		$tree->remove_unbranched_internals;
	}

	my $ntax_pruned = $tree->get_ntax;
	$logger->info("Taxon number changed from $ntax_initial to $ntax_pruned");
	
	# write rerooted tree to output file
	open my $out, '>', $outfile or die $!;	
	print $out $tree->to_newick(nodelabels=>1);
	close $out;
	
	$logger->info("DONE. Outfile written to $outfile");
	return 1;
}

1;
