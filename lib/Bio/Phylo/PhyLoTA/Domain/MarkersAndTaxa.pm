# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use List::MoreUtils 'uniq';
use Bio::Phylo::Matrices::Datum;


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

=item calc_mean_distance

Calculates the average pairwise distance within the alignment.

=cut

sub calc_mean_distance {
	my ( $class, $fastastr ) = @_;
	my $stream = Bio::AlignIO->new(-string=>$fastastr, -format=>"fasta");  
	my $aln = $stream->next_aln;  
	return (1 - $aln->average_percentage_identity/100)
}

=item dedup

Removes duplicate sequences (by GI)

=cut

sub dedup {
	my ( $class, %fasta ) = @_;
	my %seen;
	my %result;
	for my $defline ( keys %fasta ) {
		if ( $defline =~ /gi\|(\d+)/ ) {
			my $gi = $1;
			$result{$defline} = $fasta{$defline} unless $seen{$gi}++;			
		}
	}
	return %result;
}

=item parse_taxa_file

Reads the provided taxa file, returns a list of records where each consists of a 
hash where keys are the headers in the taxa table.

B<Note>: this method, and C<parse_taxa_string> are completely agnostic about the
column headers so the can equally be used to read other tab-separated files, such
as the backbone markers table, for example.

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
		
		# skip comments and blank lines
		next LINE if /^\s*#/;
		next LINE if /^\s*$/;
		
		# store header
		if ( not @header ) {
			@header = split /\t/, $line;
			next LINE;
		}
		my %record = ( 'keys' => \@header );
		
		# store fields
		my @fields = split /\t/, $line;
		for my $i ( 0 .. $#header ) {
			$record{$header[$i]} = $fields[$i];
		}
		
		# store record
		push @result, \%record;
	}
	return @result;
}

=item get_supermatrix_numtaxa

Given the file location of a supermatrix, returns the number of taxa

=cut

sub get_supermatrix_numtaxa {
	my ($self, $supermatrix) = @_;
	# parse number of sites in alignment from first line in supermatrix
	open my $file, '<', $supermatrix or die $!;
	my $firstline = <$file>;
	$firstline =~ s/^\s+//;
	my $numsites = (split(/\s/, $firstline))[0];
	close $file;
	return $numsites
}

=item get_supermatrix_numsites

Given the file location of a supermatrix, returns the number of sites

=cut

sub get_supermatrix_numsites {
	my ($self, $supermatrix) = @_;
	# parse number of sites in alignment from first line in supermatrix
	open my $file, '<', $supermatrix or die $!;
	my $firstline = <$file>;
	$firstline =~ s/^\s+//;
	my $numsites = (split(/\s/, $firstline))[1];
	close $file;
	return $numsites
}


=item get_distinct_taxa

Returns the distinct IDs for the provided taxonomic level.

=cut

sub get_distinct_taxa {
	my ( $class, $taxon, @records ) = @_;
	my %ids = map { $_ => 1 } map { $_->{$taxon} } @records;
	# filter for e.g. 'NA' values
	my @taxa = grep /[0-9]+/,  keys %ids;
	return @taxa;
}

=item get_highest_informative_level

Determines, for a given species table, which taxonomic level
(column in the table) is the highest one that can give information about
distinct taxa in this table. This will be the column which does not contain NA values and 
which has at least two different taxon identifiers. This level forms a good basis to check
for paraphyly in the branches spawning from there.

=cut

sub get_highest_informative_level{
	my ( $class, @records) = @_;	
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $result;
	
	# get all possible taxonomic ranks, ordered from highest to lowest
	my @all_ranks = $mts->get_taxonomic_ranks;
	my %ranks_in_table = map { $_=>1 } keys %{$records[0]};
	
	# iterate over all ranks; the first informative taxon level is the level that 
	# does not contain NA values and that has distinct taxon IDs
	foreach my $rank (@all_ranks){		

		# skip rank if not in our taxa table
		next if not exists $ranks_in_table{$rank};
		
		# extract column for taxonomic level
		my @column = map {my %h=%$_; my $entry=$h{$rank}; $entry=~s/\s//g; $entry; } @records;
	
		# omit columns that contain 'NA' values
		next if grep {$_ eq 'NA'} @column;
	
		if (scalar (uniq @column) > 1){
			$result = $rank;
			last;
		}
	}
	return $result;
}

=item get_root_taxon_level

Returns the level of the 'root' taxon. This is the lowest taxonomic rank which
is the same for all entries in the species table.

=cut

sub get_root_taxon_level {
	my ( $class, @records) = @_;	
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $result;
	
	# get all possible taxonomic ranks, ordered from lowest to highest
	my @all_ranks = reverse $mts->get_taxonomic_ranks;
	my %ranks_in_table = map { $_=>1 } keys %{$records[0]};
	
	# iterate over all ranks; the root taxon level is the first rank for which all 
	#  entries in the taxa table are the same!
	foreach my $rank (@all_ranks){		

		# skip rank if not in our taxa table
		next if not exists $ranks_in_table{$rank};
		
		# extract column for taxonomic level
		my @column = map {my %h=%$_; my $entry=$h{$rank}; $entry=~s/\s//g; $entry; } @records;
	
		# omit columns that contain 'NA' values
		next if grep {$_ eq 'NA'} @column;
		if (scalar (uniq @column) == 1){
			$result = $rank;
			last;
		}
	}
	return $result;
}

=item query_taxa_table

given a taxon id or a reference to a list of taxon identifiers, and a single taxonomic rank or a reference to a list of ranks,
returns all entries in the taxon table matching the rank(s) that are associated
with the input taxon ID(s). 

=cut

sub query_taxa_table {
	my ($self, $tids, $ranks, @records) = @_;
	
	# input can either be scalars or array refs
	my %tids = ref($tids) ? map {$_=>1} @${tids} : ($tids=>1);
	my %ranks = ref($ranks) ? map {$_=>1} @${ranks} : ($ranks=>1);
	
	my @matching = grep { my %h = %{$_}; grep{ exists($tids{$_}) } values%h } @records;	
	my @ids = uniq map { @{$_}{keys(%ranks)} } @matching;
	# remove NA values from result
	my @result = grep /[0-9]+/, @ids;
	return @result;
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

=item get_gis_from_fasta

Extracts sequence identifiers from a hash such as is produced by parse_fasta_string.
Returns a hash ( taxon => [ gi1, gi2 ] )

=cut

sub get_gis_from_fasta {
    my ( $class, %fasta ) = @_;
    my %gis;
    for my $key ( keys %fasta ) {
        if ( $key =~ /taxon\|(\d+)/ ) {
            my $taxon = $1;
            $gis{$taxon} = [] if not $gis{$taxon};
            if ( $key =~ /(>|\|)gi\|(\d+)/ ) {
            	my $gi = $1;
            	push @{ $gis{$taxon} }, $gi;
            }
        }
    }
    return %gis;

}

=item get_sequences_for_taxon

Extracts a list of sequence strings for the provided taxon from a hash such as is
produced by parse_fasta_string

=cut

sub get_sequences_for_taxon {
	my ( $class, $taxon, %fasta ) = @_;
	my @sequences;
	for my $key ( sort keys %fasta ) {
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