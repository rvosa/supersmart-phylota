# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Service::CalibrationTableCreator;
use strict;
use warnings;
use Moose;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Domain::CalibrationTable;

extends 'Bio::Phylo::PhyLoTA::Service';

=head1 NAME

Bio::Phylo::PhyLoTA::Service::CalibrationTableCreator - creates sets of calibration points

=head1 DESCRIPTION

Creates a L<Bio::Phylo::PhyLoTA::Domain::CalibrationTable> object by placing
L<Bio::Phylo::PhyLoTA::Domain::FossilData> objects on a taxonomy.

=head1 METHODS

=over

=item new

The constructor takes no arguments.

=cut

sub new {
    my $class = shift;
    my $self = bless {}, $class;
    return $self;
}

=item find_calibration_point

Given a L<Bio::Phylo::PhyLoTA::Domain::FossilData>, queries for any nested instances
of genus and family in the taxonomy that match those of the fossil. If found, attaches
the genus as the point which the fossil calibrates.

=cut

sub find_calibration_point {
    my ( $self, $fd ) = @_;
    my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
    
    # search for the family and the genus of the fossil datum
    my ( $genus, $family );    
    ( $family ) = $mts->get_nodes_for_names( $fd->family ) if $fd->family;
    ( $genus  ) = $mts->get_nodes_for_names( $fd->genus  ) if $fd->genus;
    
    # if they are indeed nested and it's a stem fossil, attach the genus
    # node as the calibration point
    if ( $family && $genus && $family->is_ancestor_of($genus) && $fd->crown_stem eq 'stem' ) {
        $fd->point( $genus );
    }
    return $fd;
}

=item create_calibration_table

Given a tree and an array of fossils, creates a 
L<Bio::Phylo::PhyLoTA::Domain::CalibrationTable>.

=cut

sub create_calibration_table {
    my ( $self, $tree, @fossildata ) = @_;
    my $table = Bio::Phylo::PhyLoTA::Domain::CalibrationTable->new;
    for my $fd ( @fossildata ) {
        if ( my $node = $fd->point ) {
            
            # fetch all the tips in the taxonomy using a level-order,
            # non-recursive, queue-based traversal
            my @tips;
            my @queue = ($node);
            while ( @queue ) {
                my $node = shift @queue;
                my @children = @{ $node->get_children };
                if ( @children ) {
                    push @queue, @children;
                }
                else {
                    
                    # search for the equivalent node in the tree
                    if ( my $tip = $tree->get_by_name($node->get_name) ) {
                        push @tips, $tip;
                    }
                }
            }
            
            # get most recent common ancestor
            my $mrca = $tree->get_mrca(\@tips);

            # create record in table
            $table->add_row(
                'taxa'    => [ map { $_->get_name } @{ $mrca->get_terminals } ],
                'min_age' => $fd->min_age,
                'max_age' => $fd->max_age,
            );
        }
    }
    return $table;
}

=back

=cut

1;