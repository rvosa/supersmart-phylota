#!/usr/bin/perl
use strict;
use warnings;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;

my $fasta = shift;
my %fasta   = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->parse_fasta_file($fasta);
my $ntax    = scalar keys %fasta;
my ($nchar) = map { length($_) } values %fasta;

print "$ntax $nchar\n";
for my $defline ( keys %fasta ) {
	if ( $defline =~ /taxon\|(\d+)/ ) {
		my $name = $1;
		print $name, ' ' x ( 10 - length($name) ), $fasta{$defline}, "\n";
	}
}