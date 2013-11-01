package Bio::Phylo::PhyLoTA::Domain::CalibrationPoint;
use strict;
use warnings;

our $AUTOLOAD;

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