# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Domain::FossilData;
use strict;
use warnings;

our $AUTOLOAD;

=head1 NAME

Bio::Phylo::PhyLoTA::Domain::FossilData - Fossil Data

=head1 DESCRIPTION

Object that represents a fossil datum that is instantiated from a row in a file
such as the tab-separated examples fossils1.tsv and fossils2.tsv.

=head1 METHODS

It is likely that each fossil will have at least the minimal methods that allow it
to be placed on a taxonomy, i.e. the C<calibrated_taxon>, and whether the fossil
is a crown or a stem fossil (C<crown_vs_stem>), and the age range for the fossil, i.e.
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

=item fossil_name

Returns the name of the fossil from the fossil table.
If provided, fossil_name it is set to the given value.

=cut

sub fossil_name {
    my $self = shift;
    $self->{'FossilName'} = shift if @_;
    return $self->{'FossilName'};
}

=item crown_vs_stem

Returns either 'crown' or 'stem' depending 
on the fossil entry in the provided fossil table.  
If provided, 'CrownvsStem' is set to the given value for this object.

=cut

sub crown_vs_stem {
    my $self = shift;
    $self->{'CrownvsStem'} = shift if @_;
    return $self->{'CrownvsStem'};
}

=item calibrated_taxon

Returns the name of the taxon which is calibrated according to
the fossil record.
If provided, calibrated_taxon it is set to the given value.

=cut

sub calibrated_taxon {
    my $self = shift;
    $self->{'CalibratedTaxon'} = shift if @_;
    return $self->{'CalibratedTaxon'};
}

=item min_age

Getter/Setter for the minimum age of the fossil.

=cut

sub min_age {
    my $self = shift;
    $self->{'MinAge'} = shift if @_;
    return $self->{'MinAge'};
}

=item best_practice_score

Getter/Setter for the fossil best practice score.

=cut

sub best_practice_score {
    my $self = shift;
    $self->{'BestPracticeScore'} = shift if @_;
    return $self->{'BestPracticeScore'};
}

=item fossil_author

Getter/Setter for the author of the record (NUsrcrFos column in the table).

=cut

sub fossil_author {
    my $self = shift;
    $self->{'NUsrcrFos'} = shift if @_;
    return $self->{'NUsrcrFos'};
}

=item date_added

Getter/Setter for the date the fossil record was added.

=cut

sub date_added {
    my $self = shift;
    $self->{'DcrFos'} = shift if @_;
    return $self->{'DcrFos'};
}

=item calibration_points

Returns an array with the calibration points found for the fossil record
according to the taxon name of the record. Note that there can be more than one 
possible calibration point, since taxa name are not unique. 

=cut

sub calibration_points {
    my $self = shift;    
    $self->{'calibration_points'} = [@_] if @_;
    return defined($self->{'calibration_points'}) ? @{$self->{'calibration_points'}} : ();
}

1;
