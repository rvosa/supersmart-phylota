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
'min_age', 'max_age' and 'taxa' (i.e. the names of the subtended terminal taxa). These
properties subsequently become available as object methods, i.e. $point->min_age and
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
			
1;