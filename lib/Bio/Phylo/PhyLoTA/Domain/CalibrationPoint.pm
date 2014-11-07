package Bio::Phylo::PhyLoTA::Domain::CalibrationPoint;
use strict;
use warnings;

our $AUTOLOAD;

=head1 NAME

Bio::Phylo::PhyLoTA::Domain::CalibrationPoint - a single, known node age

=head1 DESCRIPTION

Represents a single record in a set of calibration points.

=head1 METHODS

=over

=item new

The constructor is typically executed by the 
L<Bio::Phylo::PhyLoTA::Domain::CalibrationTable>, which passes it the named arguments
'min_age' and/or 'max_age' and 'taxa' (i.e. the names of the subtended terminal taxa). 
These properties subsequently become available as object methods, i.e. $point->min_age and
so on.

=back

=cut

sub new {
	my $class = shift;
	my %args = @_; 
	my $self = \%args;
	return bless $self, $class;
}
	
sub AUTOLOAD {
	my $self = shift;
	my $method = $AUTOLOAD;
	$method =~ s/.+://;
	if ( $method !~ /^[A-Z]$/ ) {
		my $value = $self->{$method};
		if ( ref $value and ref $value eq 'ARRAY' ) {
			return @{ $value };
		} 
		else {
			return $value;
		}
	}
}

=item min_age

    Getter/setter for minimum age for a L<Bio::Phylo::PhyLoTA::Domain::CalibrationPoint> object.
    Note that this corresponds to the column named "MinAge" in the fossil table.

=cut 

sub min_age {
        my $self = shift;
        if ( @_ ) {
                $self->{'min_age'} = shift;
        }
        return $self->{'min_age'};
}

=item max_age

    Getter/setter for maximum age for a L<Bio::Phylo::PhyLoTA::Domain::CalibrationPoint> object.
    Note that this corresponds to the column named "MaxAge" in the fossil table.

=cut 

sub max_age {
        my $self = shift;
        if ( @_ ) {
                $self->{'max_age'} = shift;
        }
        return $self->{'max_age'};
}

=item name

    Getter/setter for the calibrated taxon name  for a  L<Bio::Phylo::PhyLoTA::Domain::CalibrationPoint> object.
    Note that this corresponds to the column named "CalibratedTaxon" in the fossil table.

=cut 

sub name {
        my $self = shift;
        if ( @_ ) {
                $self->{'name'} = shift;
        }
        return $self->{'name'};
}


			
1;
