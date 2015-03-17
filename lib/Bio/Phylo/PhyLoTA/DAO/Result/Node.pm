# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::DAO::Result::Node;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

Bio::Phylo::PhyLoTA::DAO::Result::Node

=cut

__PACKAGE__->table("nodes");

=head1 ACCESSORS

=head2 ti

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 ti_anc

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 terminal_flag

  data_type: 'tinyint'
  is_nullable: 1

=head2 rank_flag

  data_type: 'tinyint'
  is_nullable: 1

=head2 model

  data_type: 'tinyint'
  is_nullable: 1

=head2 taxon_name

  data_type: 'varchar'
  is_nullable: 1
  size: 128

=head2 common_name

  data_type: 'varchar'
  is_nullable: 1
  size: 128

=head2 rank

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 n_gi_node

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_gi_sub_nonmodel

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_gi_sub_model

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_clust_node

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_clust_sub

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_piclust_sub

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_sp_desc

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_sp_model

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_leaf_desc

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_otu_desc

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "ti",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "ti_anc",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "terminal_flag",
  { data_type => "tinyint", is_nullable => 1 },
  "rank_flag",
  { data_type => "tinyint", is_nullable => 1 },
  "model",
  { data_type => "tinyint", is_nullable => 1 },
  "taxon_name",
  { data_type => "varchar", is_nullable => 1, size => 128 },
  "common_name",
  { data_type => "varchar", is_nullable => 1, size => 128 },
  "rank",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "n_gi_node",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_gi_sub_nonmodel",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_gi_sub_model",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_clust_node",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_clust_sub",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_piclust_sub",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_sp_desc",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_sp_model",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_leaf_desc",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_otu_desc",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);
__PACKAGE__->set_primary_key("ti");

=head1 RELATIONS

=head2 clusters

Type: has_many

Related object: L<Bio::Phylo::PhyLoTA::DAO::Result::Cluster>

=cut

__PACKAGE__->has_many(
  "clusters",
  "Bio::Phylo::PhyLoTA::DAO::Result::Cluster",
  { "foreign.ti_root" => "self.ti" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 specimens

Type: has_many

Related object: L<Bio::Phylo::PhyLoTA::DAO::Result::Specimen>

=cut

__PACKAGE__->has_many(
  "specimens",
  "Bio::Phylo::PhyLoTA::DAO::Result::Specimen",
  { "foreign.ti" => "self.ti" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07002 @ 2012-05-29 21:38:33
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:gASJ9/nj5PXGPX9jZT2prA


# You can replace this text with custom content, and it will be preserved on regeneration
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Forest::NodeRole;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
push @Bio::Phylo::PhyLoTA::DAO::Result::Node::ISA, 'Bio::Phylo::Forest::NodeRole';

my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
my %tree;

=head2 table

Getter/setter that maps this ORM class onto the correct version (184) of the underlying
database table.

=cut

sub table {
	my $class = shift;
	my $table = shift;
	my $release = Bio::Phylo::PhyLoTA::Config->new->currentGBRelease;
	$class->SUPER::table( $table . '_' . $release );
}

=head2 get_parent

Returns direct parent node.

=cut

sub get_parent {
	my $self = shift;
	if ( $self->get_generic('root') ) {
		return;
	}
	if ( my $parent_ti = $self->ti_anc ) {
		return $mts->find_node($parent_ti);
	}
	return;
}

=head2 set_parent

This is a no-op: the tree structure is immutable.

=cut

sub set_parent { return shift }

=head2 get_children

Returns array ref of direct children.

=cut

sub get_children {
	my $self = shift;
	my $ti = $self->ti;
	my @children = $mts->search_node( { ti_anc => $ti } )->all;
	return \@children;
}

=head2 get_terminal_ids

Returns array ref to all terminal children

=cut

sub get_terminal_ids {
	my $self = shift;
	my %seen;
	my @queue = @{ $self->get_children };
	while(@queue) {
		my $child = shift @queue;
		$seen{$child->ti}++;
		my @children = @{ $child->get_children };
		push @queue, @children if @children;
	}
	my @terminals = keys(%seen);
	return \@terminals;
}

=head2 get_branch_length 

Returns nothing: in this implementation (i.e. a taxonomy) there are no branch lengths

=cut

sub get_branch_length { return }

=head2 set_branch_length

This is a no-op: the tree structure is immutable.

=cut

sub set_branch_length { return shift }

=head2 get_id

Alias for C<ti>.

=cut

sub get_id { shift->ti }

=head2 set_tree

Stores a reference to the containing tree, if any.

=cut

sub set_tree {
	my ( $self, $tree ) = @_;
	$tree{ $self->get_id } = $tree;
	return $self;
}

=head2 get_tree

Returns a reference to the containing tree, if any.

=cut

sub get_tree { $tree{ shift->get_id } }

=head2 get_name

Alias for taxon_name

=cut

sub get_name { shift->taxon_name }

1;
