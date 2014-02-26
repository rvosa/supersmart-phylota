#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Temp 'tempfile';
use List::Util qw(min max sum);
use Bio::Phylo::Matrices::Datum;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;

=head1 NAME

pick_exemplars.pl - selects per-genus exemplar taxa from a set of alignments

=head1 SYNOPSIS

 $ pick_exemplars.pl --list=<list> --taxa=<taxa> [--verbose] > <supermatrix>

=head1 DESCRIPTION

Given an input file that lists alignments (one on each line), traverses each genus in
each alignment and picks the overall most divergent two species to represent their
genus. The rationale is that these species will (likely) cross the root of their genus, 
so that the genus-level tree can then be scaled to the same depth of that split and 
be grafted onto the tree without (too many) negative branch lengths.

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
my $config = Bio::Phylo::PhyLoTA::Config->new;
my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
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
		push @alignments, $_ if /\S/ && -e $_;
	}
}

# read taxa table
$log->info("going to read taxa table $taxa");
my @records = $mt->parse_taxa_file($taxa);

# extract and iterate over the distinct genera
my %genera;
for my $genus ( $mt->get_distinct_taxa( 'genus' => @records ) ) {

	# extract the distinct species for the focal genus
	my @species = $mt->get_species_for_taxon( 'genus' => $genus, @records );
	$genera{$genus} = \@species;
}
$log->info("read ".scalar(keys(%genera))." genera");

# here comes the convoluted logic:
# 1. for each alignment, calculate all pairwise distances within each genus
# 2. if multiple sequences for the same species, calculate the average distance
# 3. sort the pairs by distance within each alignment
# 4. assign the most distal pair a score that is proportional to the size of the sample
my %pairs;
my $dat = 'Bio::Phylo::Matrices::Datum';
ALN: for my $aln ( @alignments ) {
	$log->debug("assessing exemplars in $aln");
	my %fasta = $mt->parse_fasta_file($aln);
	GENUS: for my $genus ( keys %genera ) {
		next GENUS if @{ $genera{$genus} } < 2;
		$pairs{$genus} = {} if not $pairs{$genus};
	
		# lump instantiated sequence objects by species within this genus
		my %seq;
		my $sequence_count = 0;
		for my $species ( @{ $genera{$genus} } ) {
			if ( my @raw = $mt->get_sequences_for_taxon( $species, %fasta ) ) {
				$log->debug("$aln contained ".scalar(@raw)." sequences for species $species");
				my @seq = map { $dat->new( '-type' => 'dna', '-char' => $_ ) } @raw;
				$seq{$species} = \@seq;
				$sequence_count++;
			}
			else {
				$log->debug("$aln has no sequences for species $species");
			}
		}
		if ( $sequence_count == 0 ) {
			$log->debug("$aln has no sequences for genus $genus");
			next GENUS;
		}
		
		# calculate the distances within the genus
		my %dist;
		my @species = keys %seq;
		if ( scalar(@species) < 2 ) {
			$log->debug("not enough species for genus $genus in $aln");
			next GENUS;
		}
		for my $i ( 0 .. ( $#species - 1 ) ) {
			for my $j ( ( $i + 1 ) .. $#species ) {
				my $sp1 = $species[$i];
				my $sp2 = $species[$j];
				$log->debug("going to compute average distance between $sp1 and $sp2");
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
		
		# pick the most distal pair, increment its occurrence
		my ($farthest) = sort { $dist{$b} <=> $dist{$a} } keys %dist;
		$pairs{$genus}->{$farthest} += sum( 0 .. (scalar(keys(%dist))-1) );
		$log->debug("most distal pair for genus $genus in $aln: $farthest");
	}
}

# make the final set of exemplars
my @exemplars;
for my $genus ( keys %genera ) {

	# there were 2 or fewer species in the genus
	if ( not $pairs{$genus} ) {
		push @exemplars, @{ $genera{$genus} };
	}
	else {
		my %p = %{ $pairs{$genus} };
		my ($sp1,$sp2) = map { split /\|/, $_ } sort { $p{$b} <=> $p{$a} } keys %p;
		push @exemplars, $sp1, $sp2;
		$log->debug(Dumper({ $genus => \%p }));
	}
}
$log->info("identified ".scalar(@exemplars)." exemplars");

# make the best set of alignments:
# 1. map exemplar taxa to alignments and vice versa
my ( %taxa_for_aln, %alns_for_taxon, %nchar );
for my $aln ( @alignments ) {
	$taxa_for_aln{$aln} = [] if not $taxa_for_aln{$aln};
	my %fasta = $mt->parse_fasta_file($aln);
	$nchar{$aln} = $mt->get_nchar(%fasta);
	for my $taxon ( @exemplars ) {
		$alns_for_taxon{$taxon} = [] if not $alns_for_taxon{$taxon};	
		my ($seq) = $mt->get_sequences_for_taxon($taxon,%fasta);
		if ( $seq ) {
			push @{ $taxa_for_aln{$aln} }, $taxon;
			push @{ $alns_for_taxon{$taxon} }, $aln;
		}
	}
}
# 2. for each taxon, sort its alignments by decreasing taxon coverage
for my $taxon ( @exemplars ) {
	my @sorted = sort { scalar(@{$taxa_for_aln{$b}}) <=> scalar(@{$taxa_for_aln{$a}}) } @{ $alns_for_taxon{$taxon} };
	$alns_for_taxon{$taxon} = \@sorted;
}
# 3. sort the taxa by increasing occurrence in alignments
my @sorted_exemplars = sort { scalar(@{$alns_for_taxon{$a}}) <=> scalar(@{$alns_for_taxon{$b}}) } @exemplars;
my %aln;
my %seen;
TAXON: for my $taxon ( @sorted_exemplars ) {
	if ( $alns_for_taxon{$taxon} ) {
		my $aln = shift @{ $alns_for_taxon{$taxon} };
		if ( not $aln or not -e $aln ) {
			$log->warn("no alignment available for exemplar $taxon");
			next TAXON;
		}
		$aln{$aln}++;
		for my $other_taxa ( @{ $taxa_for_aln{$aln} } ) {
			delete $alns_for_taxon{$other_taxa} if ++$seen{$other_taxa} == $config->BACKBONE_MIN_COVERAGE;
		}
	}
	last TAXON unless scalar keys %alns_for_taxon;
}
my @filtered_exemplars = grep { ! $alns_for_taxon{$_} } @sorted_exemplars;

# produce interleaved phylip output
my @sorted_alignments = keys %aln;
my $nchar = sum( @nchar{@sorted_alignments} );
my $ntax  = scalar(@filtered_exemplars);
my $names_printed;
print $ntax, ' ', $nchar, "\n";
for my $aln ( @sorted_alignments ) {
	$log->info("including $aln in supermatrix");
	my %fasta = $mt->parse_fasta_file($aln);
	for my $taxon ( @filtered_exemplars ) {
		if ( not $names_printed ) {
			print $taxon, ' ' x ( 10 - length($taxon) );
		}
		my ( $seq ) = $mt->get_sequences_for_taxon($taxon,%fasta);
		if ( $seq ) {
			print $seq;
		}
		else {
			print '?' x $nchar{$aln};
		}
		print "\n";
	}
	$names_printed++;
}
