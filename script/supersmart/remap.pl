#!/usr/bin/perl
use strict;
use warnings;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Service;

my $infile = shift;
my $service = Bio::Phylo::PhyLoTA::Service->new;

my $tree = parse_tree(
	'-format'     => 'newick',
	'-file'       => $infile,
	'-as_project' => 1,
);

$tree->visit(sub{
	my $n = shift;
	if ( $n->is_terminal ) {
		my $id = $n->get_name;
		$n->set_name( $service->find_node($id)->taxon_name );
	}
});

print $tree->to_newick;