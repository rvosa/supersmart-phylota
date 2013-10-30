#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;

# process command line arguments
my $verbosity = WARN;
my ( @taxa, $alignments, $locus, $suffix );
GetOptions(
	'taxa=i'       => \@taxa,
	'alignments=s' => \$alignments,
	'locus=s'      => \$locus,
	'suffix=s'     => \$suffix,
	'verbose+'     => \$verbosity,
);

# instantiate helper objects
my $mts = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
my $sg  = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# fetch candidate files from workdir
my @files;
{
	$log->info("going to fetch candidate alignments from $alignments");
	open my $fh, '<', $alignments or die $!;
	@files = grep { /\.\Q$locus\E\./ } grep { /\.fa$/ } <$fh>; 
	chomp @files;
	close $fh;
	$log->info("found ".scalar(@files)." candidates in $alignments");
}

# fetch occurrence counts for each file
my %counts;
for my $file ( @files ) {
	my %taxa = map { $_ => 1 } @taxa;
	my $count;
	open my $fh, '<', $file or die "Can't open ${file}: $!";
	while(<$fh>) {
		next unless /^>/;
		if ( /taxon\|(\d+)/ ) {
			my $taxon = $1;
			$count++ if $taxa{$taxon};
		}
	}
	$counts{$file} = $count || 0;
	$log->info("file $file has $counts{$file} occurrences of the taxa of interest");
}

# sort files
my @sorted = sort { $counts{$b} <=> $counts{$a} } @files;

# write the ones of interest
FILE: for my $file ( @sorted ) {
	if ( $file =~ /(\d+)\.[^\/]*\.fa$/ ) {
		my $gi = $1;
		my %fasta   = $mts->parse_fasta_file( $file );
		my ( $seq ) = map { $fasta{$_} } grep { /gi\|$gi\|/ } keys %fasta;
		my @indices = $sg->get_aligned_locus_indices( $gi, $locus, $seq );

		# extract subsequences
		for my $defline ( keys %fasta ) {
			my @seq = split //, $fasta{$defline};
			my @subseq = @seq[@indices];
			$fasta{$defline} = join '', @subseq;
		}

		# now keep taxa of interest
		my %reduced = $mts->keep_taxa( \@taxa, \%fasta );
		%reduced = degap(%reduced);

		# skip once we are below 3 taxa of interest
		if ( 3 <= scalar keys %reduced ) {
			my $outfile = $file;
			$outfile =~ s/\.fa$/.$suffix.fa/;
			open my $fh, '>', $outfile or die $!;
			for my $defline ( keys %reduced ) {
				print $fh '>', $defline, "\n";
				print $fh $reduced{$defline}, "\n";
			}
		}
		else {
			$log->info("file $file has too few sequences of interest, skipping");
			next FILE;
		}
	}
}

sub degap {
	my %fasta  = @_;
	my @matrix = map { [ split //, $_ ] } values %fasta;
	my $nchar = $#{ $matrix[0] };
	my %indices;
	SITE: for my $i ( 0 .. $nchar ) {
		for my $row ( @matrix ) {
			$indices{$i} = 1;
			next SITE if $row->[$i] ne '-';
		}
		delete $indices{$i};
	} 
	my @keep = sort { $a <=> $b } keys %indices;
	for my $defline ( keys %fasta ) {
		my @seq = split //, $fasta{$defline};
		my @subseq = @seq[@keep];
		$fasta{$defline} = join '', @subseq;
	}
	return %fasta;
}
