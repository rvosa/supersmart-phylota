# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Domain::Sequence;
use strict;
use warnings;

=head1 NAME

Bio::Phylo::PhyLoTA::Domain::Sequence

=head1 DESCRIPTION

Table with sequences. XXX: kill me.

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

