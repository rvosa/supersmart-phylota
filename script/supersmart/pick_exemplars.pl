#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);
use Bio::Phylo::Matrices::Datum;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;

=head1 NAME

pick_exemplars.pl - selects per-genus exemplar taxa from a set of alignments

=head1 SYNOPSIS

 $ pick_exemplars.pl --list=<list> --taxa=<taxa> [--verbose] > <supermatrix>

=head1 DESCRIPTION

Given an input file that lists alignments (one on each line), traverses each genus in
each alignment and picks the overall most divergent two species to represent their
genus. The rationale is that these species will cross the root of their genus, so that
the genus-level tree can then be scaled to the same depth and grafted onto the tree.

B<Point of consideration>: the node depths on the exemplar tree will be underestimates
relative to the genus-level tree (due to the node density effect), so it might be better
to give the node depths from the exemplar tree as priors and then see where the posteriors
come out. Disadvantage is that this is likely to lead to negative branch lengths in some
cases.

=cut

# process command line arguments
my $verbosity = WARN;
my ( $list, $taxa );
GetOptions(
	'list=s'   => \$list,
	'taxa=s'   => \$taxa,
	'verbose+' => \$verbosity,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-class' => 'main',
	'-level' => $verbosity,
);

# read list of alignments
my @alignments;
{
	$log->info("going to read alignment list $list");
	open my $fh, '<', $list or die $!;
	while(<$fh>) {
		chomp;
		push @alignments, $_;
	}
}

# read taxa table
my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
my @records = $mt->parse_taxa_file($taxa);

# extract and iterate over the distinct genera
my %genera;
for my $genus ( $mt->get_distinct_taxa( 'genus' => @records ) ) {

	# extract the distinct species for the focal genus
	my @species = $mt->get_species_for_taxon( 'genus' => $genus, @records );
	$genera{$genus} = \@species;
}

# here comes the convoluted logic:
# 1. for each alignment, calculate all pairwise distances within each genus
# 2. if multiple sequences for the same species, calculate the average distance
# 3. normalize the distances, e.g. between 0 and 1. This is so that we don't pick
#    as exemplar pair two taxa that just happened to be divergent in a rapidly
#    evolving locus even though they are not the overall most divergent
# 4. sort the distances, retaining the most distal pair for each alignment
# 5. pick the modal pair across alignments
my %pairs;
my $dat = 'Bio::Phylo::Matrices::Datum';
for my $aln ( @alignments ) {
	my %fasta = $mt->parse_fasta_file($aln);
	GENUS: for my $genus ( keys %genera ) {
		next GENUS if @{ $genera{$genus} } < 3;	
		$pairs{$genus} = {} if not $pairs{$genus};
	
		# lump instantiated sequence objects by species within this genus
		my %seq;
		for my $species ( @{ $genera{$genus} } ) {
			if ( my @raw = $mt->get_sequences_for_taxon( $species, %fasta ) ) {
				$log->info("$aln contained ".scalar(@raw)." sequences for $species");
				my @seq = map { $dat->new( '-type' => 'dna', '-char' => $_ ) } @raw;
				$seq{$species} = \@seq;
			}
			$log->info("$aln has no sequences for $species");
		}
		
		# calculate the distances within the genus
		my %dist;
		my @species = keys %seq;
		for my $i ( 0 .. ( $#species - 1 ) ) {
			for my $j ( ( $i + 1 ) .. $#species ) {
				my $sp1 = $species[$i];
				my $sp2 = $species[$j];
				$log->info("going to compute average distance between $sp1 and $sp2");
				my @seqs1 = @{ $seq{$sp1} };
				my @seqs2 = @{ $seq{$sp2} };
				my $dist;
				for my $seq1 ( @seqs1 ) {
					for my $seq2 ( @seqs2 ) {
						$dist += $seq1->calc_distance($seq2);
					}
				}
				$dist /= ( scalar(@seqs1) * scalar(@seqs2) );
				my $key = join '|', sort { $a <=> $b } ( $sp1, $sp2 );
				$dist{$key} = $dist;
			}
		}
		
		# normalize the distances
		my $min = min values %dist;
		my $max = max values %dist;
		my $scale = 1 / ( $max - $min ); # XXX can divide by 0?
		for my $pair ( keys %dist ) {
			$dist{$pair} = ( $dist{$pair} - $min ) * $scale;
		}
		
		# pick the most distal pair, increment its occurrence
		my ($farthest) = sort { $dist{$b} <=> $dist{$a} } keys %dist;
		$pairs{$genus}->{$farthest}++;
	}
}

# make the final set of exemplars
my @exemplars;
for my $genus ( keys %pairs ) {

	# there were 2 or fewer species in the genus
	if ( not $pairs{$genus} ) {
		push @exemplars, @{ $genera{$genus} };
	}
	else {
		my %p = %{ $pairs{$genus} };
		my ($sp1,$sp2) = map { split /\|/, $_ } sort { $p{$b} <=> $p{$a} } keys %p;
		push @exemplars, $sp1, $sp2;
	}
}

# make the output
my %supermatrix;
for my $taxon ( @exemplars ) {
	$supermatrix{$taxon} = '';
	for my $aln ( @alignments ) {
		my %fasta = $mt->parse_fasta_file($aln);
		my $nchar = $mt->get_nchar(%fasta);
		my ($seq) = $mt->get_sequences_for_taxon($taxon,%fasta);
		if ( $seq ) {
			$supermatrix{$taxon} .= $seq;
		}
		else {
			$supermatrix{$taxon} .= '?' x $nchar;
		}
	}
}

# print the output
for my $taxon ( keys %supermatrix ) {
	if ( $supermatrix{$taxon} !~ /^\?+$/ ) {
		print '>', $taxon, "\n", $supermatrix{$taxon}, "\n";
	}
}