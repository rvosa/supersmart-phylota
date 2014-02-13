package Bio::Phylo::PhyLoTA::DAO::Result::Inparanoid;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

Bio::Phylo::PhyLoTA::DAO::Result::Inparanoid

=cut

__PACKAGE__->table("inparanoid");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 guid

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 filename

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 confidence

  data_type: 'float'
  is_nullable: 1

=head2 protid

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 bootstrap

  data_type: 'integer'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "guid",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "filename",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "confidence",
  { data_type => "float", is_nullable => 1 },
  "protid",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "bootstrap",
  { data_type => "integer", is_nullable => 1 },
);


# Created by DBIx::Class::Schema::Loader v0.07002 @ 2012-10-09 23:24:49
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:RhWtEnYbQO7qHsFdLVg9Ew


# You can replace this text with custom content, and it will be preserved on regeneration

=head2 is_orthologous

Given another Inparanoid database record, determine if the two are orthologous. This is
done by comparing the identity of their id, guid and bootstrap fields.

=cut

sub is_orthologous {
	my ( $self, $other ) = @_;
	return
		$self->id   == $other->id &&
		$self->guid == $other->guid &&
		$self->bootstrap  eq '100%' &&
		$other->bootstrap eq '100%';
}

1;
