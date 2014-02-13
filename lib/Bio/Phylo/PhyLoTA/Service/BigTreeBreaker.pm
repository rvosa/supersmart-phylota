# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Service::BigTreeBreaker;
use strict;
use warnings;

=head1 NAME

Bio::Phylo::PhyLoTA::Service::BigTreeBraker - Big Tree Braker

=head1 DESCRIPTION

Breaks tree into clades based on support, taxonomy, age, number of taxa. Returns list with 
taxonID\tdating_information. Output to be used as input to MarkerAndTaxaSelection.
XXX: kill me.

=head1 METHODS

=over

=item new

The constructor takes no arguments.

=back

=cut

sub new {
    my $class = shift;
    my $self = bless {}, $class;
    return $self;
}


1;