# this is an object oriented perl module

package Bio::SUPERSMART::Domain::CalibrationTable;
use strict;
use warnings;
use Bio::SUPERSMART::Domain::CalibrationPoint;
use Bio::Phylo::Util::Logger ':levels';

my $log = Bio::Phylo::Util::Logger->new;

=head1 NAME

Bio::SUPERSMART::Domain::CalibrationTable - Calibration table

=head1 DESCRIPTION

Object that represents a list of calibration points that can be serialized
as input for trePL.

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
L<Bio::SUPERSMART::Domain::CalibrationPoint> object.

=cut

sub add_row {
    my ( $self, %args ) = @_;
    push @{ $self }, Bio::SUPERSMART::Domain::CalibrationPoint->new(%args);
}

=item get_rows

Returns the rows in the table, which are represented as 
L<Bio::SUPERSMART::Domain::CalibrationPoint> objects.

=cut
    
sub get_rows {
    my $self = shift;
    return @{ $self };
}

=item to_string

Returns a string representation suitable for input into a tree calibration program. At 
present this simply means a serialization in the syntax that treePL uses for identifying
and dating MRCAs.

=cut

sub to_string {
    my $self = shift;
    my $string = '';
    
    # group CP's by subtended taxon set
    my %rows;
    for my $row ( $self->get_rows ) {
	my $taxa = join ',', sort { $a cmp $b } $row->taxa;
	$rows{$taxa} = [] if not $rows{$taxa};
	push @{ $rows{$taxa} }, $row;
    }

    # iterate over grouped CP's
    for my $key ( keys %rows  ) {
	my @rows = @{ $rows{$key} };
        my $row;

	# pick the "right" on if multiple
	if ( scalar(@rows) > 1 ) {
            my $names = join ', ', map { $_->nfos } @rows;
	    $log->warn("Multiple calibration points for same node ($names), choosing oldest");
            ($row) = sort { $b->max_age <=> $a->max_age } @rows;
	}
	else {
	    ($row) = @rows;
	}

	# write output
        my @taxa = $row->taxa;
        $string .= "mrca = NFos" . $row->nfos . " @taxa\n";
        if ( $row->max_age ) {
            $string .= "max = NFos" . $row->nfos . " " .  $row->max_age . "\n";
        }
        if ( $row->min_age ){
            $string .= "min = NFos" . $row->nfos . " " . $row->min_age . "\n";
        }               
    }
    return $string;
}

=item sort_by_min_age

Sorts the table according to the min_age attributes of 
all L<Bio::SUPERSMART::Domain::CalibrationPoint> objects contained in the
table. Sorting is done in ascending order. This function also sets prefixes in front 
of the calibrated taxon name for each L<Bio::SUPERSMART::Domain::CalibrationPoint>, 
such that also the names have ascending order. This function is to get control
over the behavior of treePL, which seems to care about the 
order of the calibrated taxa recorded in the config file.

=cut

sub sort_by_min_age {
        my $self = shift;       

        # sort by ascending min taxon age
        @{$self} = sort { $a->min_age <=> $b->min_age } $self->get_rows;                         

        # loop over calibration points and prepend 'ct' and an index
        my $counter = 0;
        foreach my $cp ( @{$self} ) {
                my $curr_name = $cp->name;
                # ct stands for 'calibrated taxon'
                $cp->name('ct'.$counter++.'_'.$curr_name)
        }
        return;
}

=item remove_orphan_taxa

Removes calibration points from the table that only contain a single leaf taxon to be calibrated. 
Since internal nodes of SUPERSMART trees are not assigned to any taxon id, any taxon id in the
calibration table has rank species (or lower). Calibration is done by selecting the mrca of all leave taxa
for one calibration point and setting the age of that internal node to the corresponding age.
If there is only one leave for a calibration point, one would calibrate the leave representing an
extant species. This functin removes these 'orphan' calibration points from the table. 

=cut


sub remove_orphan_taxa {
        my $self = shift;
        @{$self} = grep { $_->taxa > 1 } $self->get_rows;
        return;
}

=back

=cut

1;
