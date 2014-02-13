# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Domain::CalibrationTable;
use strict;
use warnings;
use Bio::Phylo::PhyLoTA::Domain::CalibrationPoint;

=head1 NAME

Bio::Phylo::PhyLoTA::Domain::CalibrationTable - Calibration table

=head1 DESCRIPTION

Object that represents a list of calibration points that can be serialized
as input for treePL.

=head1 METHODS

=over

=item new

The constructor takes no arguments.

=cut

sub new {
    my $class = shift;
    my $self = bless [], $class;
    return $self;
}

=item add_row

Adds a row to the table. Each row is represented as a 
L<Bio::Phylo::PhyLoTA::Domain::CalibrationPoint> object.

=cut

sub add_row {
    my ( $self, %args ) = @_;
    push @{ $self }, Bio::Phylo::PhyLoTA::Domain::CalibrationPoint->new(%args);
}

=item get_rows

Returns the rows in the table, which are represented as 
L<Bio::Phylo::PhyLoTA::Domain::CalibrationPoint> objects.

=cut
    
sub get_rows {
    my $self = shift;
    return @{ $self };
}

=item to_string

Returns a string representation suitable for input into
a tree calibration program. XXX: this needs to be fixed/enhanced to be able
to provide the right representation for multiple programs (e.g. *BEAST, treePL).

=back

=cut

sub to_string {
    my $self = shift;
    my $string = '';
    for my $row ( @{ $self } ) {
        $string .= join(' ', @{ $row->{taxa} }, '|', $row->{min_age}, $row->{max_age} ) . "\n";
    }
    return $string;
}


1;