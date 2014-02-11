# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Domain::FossilData;
use strict;
use warnings;

our $AUTOLOAD;

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

=head1 NAME

Bio::Phylo::PhyLoTA::Domain::FossilData - Fossil Data

=head1 DESCRIPTION

Object that represents a fossil datum that is instantiated from a row in a file
such as in $config->FOSSIL_TABLE_FILE

=cut

