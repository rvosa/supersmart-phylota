# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use strict;
use warnings;


=head1 NAME

Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa - Markers and Taxa

=head1 DESCRIPTION

This package provides some static methods for interacting with FASTA files. XXX: this
needs a different name and should be moved into a service package.

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

=item parse_fasta_file

Reads the provided FASTA file, returns a hash where keys are FASTA definition lines,
values are sequences.

=cut

sub parse_fasta_file {
    my ( $class, $file ) = @_;
    open my $fh, '<', $file or die $!;
    my $string = do { local $/; <$fh> };
    return $class->parse_fasta_string($string);
}

=item parse_fasta_string

Reads the provided FASTA string, returns a hash where keys are FASTA definition lines,
values are sequences.

=cut

sub parse_fasta_string {
    my ( $class, $string ) = @_;
    my @lines = split /\n/, $string;
    my %fasta;
    my $current;
    for my $line ( @lines ) {
        chomp $line;
        if ( $line =~ /^>(.+)/ ) {
            $current = $1;
            if ( exists $fasta{$current} ) {
                $fasta{$current} = '';
            }            
        }
        else {
            $fasta{$current} .= $line;
        }
    }
    return %fasta;    
}

=item parse_taxa_file

Reads the provided taxa file, returns a list of records where each consists of a 
hash where keys are the headers in the taxa table.

=cut

sub parse_taxa_file {
	my ( $class, $file ) = @_;
	open my $fh, '<', $file or die $!;
	my $string = do { local $/; <$fh> };
	return $class->parse_taxa_string($string);
	
}

=item parse_taxa_string

Reads the provided taxa string, returns a list of records where each consists of a 
hash where keys are the headers in the taxa table.

=cut

sub parse_taxa_string {
	my ( $class, $string ) = @_;
	my @header;
	my @result;
	LINE: for my $line ( split /\n/, $string ) {
		chomp $line;
		if ( not @header ) {
			@header = split /\t/, $line;
			next LINE;
		}
		my @fields = split /\t/, $line;
		my %record;
		for my $i ( 0 .. $#header ) {
			$record{$header[$i]} = $fields[$i];
		}
		push @result, \%record;
	}
	return @result;
}

=item get_distinct_taxa

Returns the distinct IDs for the provided taxonomic level.

=cut

sub get_distinct_taxa {
	my ( $class, $taxon, @records ) = @_;
	my %ids = map { $_ => 1 } map { $_->{$taxon} } @records;
	return keys %ids;
}

=item get_species_for_taxon

Returns the distinct species IDs for the provided taxon ID.

=cut

sub get_species_for_taxon {
	my ( $class, $taxon, $id, @records ) = @_;
	my %ids = map { $_ => 1 } map { $_->{'species'} } grep { $_->{$taxon} == $id } @records;
	return keys %ids;
}

=item get_taxa_from_fasta

Extracts taxon identifiers from a hash such as is produced by parse_fasta_string

=cut

sub get_taxa_from_fasta {
    my ( $class, %fasta ) = @_;
    my @taxa;
    for my $key ( keys %fasta ) {
        if ( $key =~ /taxon\|(\d+)/ ) {
            my $taxon = $1;
            push @taxa, $taxon;
        }
    }
    return @taxa;
}

=item get_sequences_for_taxon

Extracts a list of sequence strings for the provided taxon from a hash such as is
produced by parse_fasta_string

=cut

sub get_sequences_for_taxon {
	my ( $class, $taxon, %fasta ) = @_;
	my @sequences;
	for my $key ( keys %fasta ) {
		if ( $key =~ /taxon\|$taxon[^\d]/ ) {
			push @sequences, $fasta{$key};
		}
	}
	return @sequences;
}

=item get_nchar

Returns the number of characters in the fasta hash.

=cut

sub get_nchar {
	my ( $class, %fasta ) = @_;
	my ($seq) = values %fasta;
	return length $seq;
}

=item keep_taxa

Retains the records for the provided taxon identifiers that occur in a hash such as is 
produced by parse_fasta_string

=cut

sub keep_taxa {
	my ( $class, $taxa, $fasta ) = @_;
	my %result;
	my %taxa = map { $_ => 1 } @{ $taxa };
	for my $defline ( keys %{ $fasta } ) {
		if ( $defline =~ /taxon\|(\d+)/ ) {
			my $taxon = $1;
			$result{$defline} = $fasta->{$defline} if $taxa{$taxon};
		}
	}
	return %result;
}

=back

=cut

1;