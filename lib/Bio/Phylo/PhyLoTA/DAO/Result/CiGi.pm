# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::DAO::Result::CiGi;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

Bio::Phylo::PhyLoTA::DAO::Result::CiGi

=cut

__PACKAGE__->table("ci_gi");

=head1 ACCESSORS

=head2 ti

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 clustid

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 cl_type

  data_type: 'enum'
  extra: {list => ["node","subtree"]}
  is_nullable: 1

=head2 gi

  data_type: 'bigint'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 ti_of_gi

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "ti",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "clustid",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "cl_type",
  {
    data_type => "enum",
    extra => { list => ["node", "subtree"] },
    is_nullable => 1,
  },
  "gi",
  { data_type => "bigint", extra => { unsigned => 1 }, is_nullable => 1 },
  "ti_of_gi",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);


# Created by DBIx::Class::Schema::Loader v0.07002 @ 2012-05-29 00:13:17
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:loJ84EFGQDeIWN/b/PoM8w


# You can replace this text with custom content, and it will be preserved on regeneration
use Bio::Phylo::PhyLoTA::Config;

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
1;
