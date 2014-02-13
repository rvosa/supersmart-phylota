# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Domain::LocalityTaxonomy;
use strict;
use warnings;

our $AUTOLOAD;

=head1 NAME

Bio::Phylo::PhyLoTA::Domain::LocalityTaxonomy - Locality Taxonomy

=head1 DESCRIPTION

Table with species distribution data from GBIF (www.gbif.org). XXX: As we will not handle
the processing of GBIF data, this class will be retired.

=head1 METHODS

=over

=item new

The constructor takes as its argument a specimen.

=back

=cut

sub new {
    my $class = shift;
    my $specimen = shift;
    my $self = bless \$specimen, $class;
    return $self;
}

sub AUTOLOAD {
    my $self = shift;
    my $method = $AUTOLOAD;
    $method =~ s/.+://;
    return $$self->$method(@_);
}

1;