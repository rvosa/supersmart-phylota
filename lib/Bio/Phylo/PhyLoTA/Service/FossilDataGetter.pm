# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Service::FossilDataGetter;
use strict;
use warnings;
use Bio::Phylo::PhyLoTA::Domain::FossilData;

=head1 TITLE

Bio::Phylo::PhyLoTA::Service::FossilDataGetter - harvester of calibration data

=head1 DESCRIPTION

This module is going to fetch estimates of the age of a split in a tree from
some database. This ought to be from TimeTree, except they don't allow spiders
to crawl them, so we're going to get it from elsewhere. There is a NESCent
working group (L<http://www.nescent.org/science/awards_summary.php?id=259>)
that develops an alternative to TT. They've published an initial report
(L<http://www.citeulike.org/user/rvosa/article/10021021>) but the db is not
ready yet. The group is led by Dan Ksepka
(L<http://www4.ncsu.edu/~dtksepka/DanKsepka/Index.html>). For now, we instead
use a simple spreadsheet with fossils, but expect this to change, as per issue:
L<https://github.com/naturalis/supersmart/issues/7>

XXX: the CalibrationTableCreator and this package should probably be combined 
into something like CalibrationService

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

=item read_fossil_table

Given a fossil table, returns an array of 
L<Bio::Phylo::PhyLoTA::Domain::FossilData> objects.

=cut

sub read_fossil_table {
    my ( $self, $file ) = @_;    
    open my $fh, '<', $file or die $!;
    my ( @records, @header );
    LINE: while(<$fh>) {
        chomp;
        my @fields = split /\t/, $_;
        if ( not @header ) {
            @header = @fields;
            next LINE;
        }
        if ( @fields ) {
            my %record;
            $record{$header[$_]} = $fields[$_] for 0 .. $#header;
            push @records, Bio::Phylo::PhyLoTA::Domain::FossilData->new(\%record);
        }
    }
    close $fh;
    return @records;
}

=back

=cut

1;