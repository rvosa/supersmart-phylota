# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::DAO::Result::Cluster;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

Bio::Phylo::PhyLoTA::DAO::Result::Cluster

=cut

__PACKAGE__->table("clusters");

=head1 ACCESSORS

=head2 ti_root

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 1

=head2 ci

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 cl_type

  data_type: 'enum'
  extra: {list => ["node","subtree"]}
  is_nullable: 1

=head2 n_gi

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_ti

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 pi

  data_type: 'tinyint'
  is_nullable: 1

=head2 minlength

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 maxlength

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 maxaligndens

  data_type: 'float'
  is_nullable: 1

=head2 ci_anc

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 seed_gi

  data_type: 'bigint'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 q

  data_type: 'float'
  is_nullable: 1

=head2 tc

  data_type: 'float'
  is_nullable: 1

=head2 clustalw_tree

  data_type: 'longtext'
  is_nullable: 1

=head2 muscle_tree

  data_type: 'longtext'
  is_nullable: 1

=head2 strict_tree

  data_type: 'longtext'
  is_nullable: 1

=head2 clustalw_res

  data_type: 'float'
  is_nullable: 1

=head2 muscle_res

  data_type: 'float'
  is_nullable: 1

=head2 strict_res

  data_type: 'float'
  is_nullable: 1

=head2 ortho

  data_type: 'tinyint'
  is_nullable: 1

=head2 n_gen

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 n_child

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "ti_root",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_foreign_key => 1,
    is_nullable => 1,
  },
  "ci",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "cl_type",
  {
    data_type => "enum",
    extra => { list => ["node", "subtree"] },
    is_nullable => 1,
  },
  "n_gi",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_ti",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "pi",
  { data_type => "tinyint", is_nullable => 1 },
  "minlength",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "maxlength",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "maxaligndens",
  { data_type => "float", is_nullable => 1 },
  "ci_anc",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "seed_gi",
  { data_type => "bigint", extra => { unsigned => 1 }, is_nullable => 1 },
  "q",
  { data_type => "float", is_nullable => 1 },
  "tc",
  { data_type => "float", is_nullable => 1 },
  "clustalw_tree",
  { data_type => "longtext", is_nullable => 1 },
  "muscle_tree",
  { data_type => "longtext", is_nullable => 1 },
  "strict_tree",
  { data_type => "longtext", is_nullable => 1 },
  "clustalw_res",
  { data_type => "float", is_nullable => 1 },
  "muscle_res",
  { data_type => "float", is_nullable => 1 },
  "strict_res",
  { data_type => "float", is_nullable => 1 },
  "ortho",
  { data_type => "tinyint", is_nullable => 1 },
  "n_gen",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "n_child",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);
__PACKAGE__->set_primary_key("ci");

=head1 RELATIONS

=head2 ti_root

Type: belongs_to

Related object: L<Bio::Phylo::PhyLoTA::DAO::Result::Node>

=cut

__PACKAGE__->belongs_to(
  "ti_root",
  "Bio::Phylo::PhyLoTA::DAO::Result::Node",
  { ti => "ti_root" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "CASCADE",
    on_update     => "CASCADE",
  },
);


# Created by DBIx::Class::Schema::Loader v0.07002 @ 2012-05-29 00:13:17
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:oTSbdQg+uqkhbQlCiDai+Q


# You can replace this text with custom content, and it will be preserved on regeneration
use Bio::SUPERSMART::Config;

=head2 table

Getter/setter that maps this ORM class onto the correct version (184) of the underlying
database table.

=cut

sub table {
	my $class = shift;
	my $table = shift;
	my $release = Bio::SUPERSMART::Config->new->currentGBRelease;
	$class->SUPER::table( $table . '_' . $release );
}

=head2 clinfo

Serializes to a unique identifier string

=cut

sub clinfo {
	my $c = shift;
	return join '-', $c->seed_gi, $c->ti_root->ti, $c->ci, $c->cl_type;
}

1;
