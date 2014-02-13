# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Domain::FossilData;
use strict;
use warnings;

our $AUTOLOAD;

=head1 NAME

Bio::Phylo::PhyLoTA::Domain::FossilData - Fossil Data

=head1 DESCRIPTION

Object that represents a fossil datum that is instantiated from a row in a file
such as the tab-separated examples fossils1.tsv and fossils2.tsv. XXX: Note that, as per 
one of the milestone issues L<https://github.com/naturalis/supersmart/issues/7> this 
design may have to change.

=head1 METHODS

It is likely that each fossil will have at least the minimal methods that allow it
to be placed on a taxonomy, i.e. the C<genus> and C<family>, and whether the fossil
is a crown or a stem fossil (C<crown_stem>), and the age range for the fossil, i.e.
C<min_age> and C<max_age>. In addition, the object will have a C<point> method whose
value is a node in the taxonomy.

=over

=item new

The constructor takes no arguments.

=back

=cut

sub new {
    my $class = shift;
    my $self = shift || {};
    return bless $self, $class;
}

sub AUTOLOAD {
    my $self = shift;
    my $method = $AUTOLOAD;
    $method =~ s/.+://;
    if ( @_ ) {
        $self->{$method} = shift;
    }
    return $self->{$method};
}


1;