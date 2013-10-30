# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use strict;
use warnings;

sub new {
    my $class = shift;
    my $self = bless {}, $class;
    return $self;
}

sub parse_fasta_file {
    my ( $class, $file ) = @_;
    open my $fh, '<', $file or die $!;
    my $string = do { local $/; <$fh> };
    return $class->parse_fasta_string($string);
}

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


1;

=head1 NAME

Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa - Markers and Taxa

=head1 DESCRIPTION

Table of markers and their taxa. 

=cut

