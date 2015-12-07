package Bio::Phylo::PhyLoTA::DAO::Result::Specimen;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

Bio::Phylo::PhyLoTA::DAO::Result::Specimen

=cut

__PACKAGE__->table("specimens");

=head1 ACCESSORS

=head2 si

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 ti

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 1

=head2 latitude

  data_type: 'float'
  is_nullable: 1

=head2 longitude

  data_type: 'float'
  is_nullable: 1

=head2 min_age_years

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 max_age_years

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "si",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "ti",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_foreign_key => 1,
    is_nullable => 1,
  },
  "latitude",
  { data_type => "float", is_nullable => 1 },
  "longitude",
  { data_type => "float", is_nullable => 1 },
  "min_age_years",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "max_age_years",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);
__PACKAGE__->set_primary_key("si");

=head1 RELATIONS

=head2 ti

Type: belongs_to

Related object: L<Bio::Phylo::PhyLoTA::DAO::Result::Node>

=cut

__PACKAGE__->belongs_to(
  "ti",
  "Bio::Phylo::PhyLoTA::DAO::Result::Node",
  { ti => "ti" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "CASCADE",
    on_update     => "CASCADE",
  },
);


# Created by DBIx::Class::Schema::Loader v0.07002 @ 2012-05-29 21:38:33
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:c2JeyXmicK57u5YJc/WmWg


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
1;
