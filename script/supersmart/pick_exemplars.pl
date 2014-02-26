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
my $dat    = 'Bio::Phylo::Matrices::Datum';
my $config = Bio::Phylo::PhyLoTA::Config->new;
my $mt     = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
my $log    = Bio::Phylo::Util::Logger->new(
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

# extract the distinct species for each genus
my %species_for_genus;
for my $genus ( $mt->get_distinct_taxa( 'genus' => @records ) ) {

	# extract the distinct species for the focal genus
	my @species = $mt->get_species_for_taxon( 'genus' => $genus, @records );
	$species_for_genus{$genus} = \@species;
}
$log->info("read ".scalar(keys(%species_for_genus))." genera");

# this will store distances within each genus
my %distance = map { $_ => {} } keys %species_for_genus;

# iterate over alignments
ALN: for my $aln ( @alignments ) {
	$log->debug("assessing exemplars in $aln");
	my %fasta = $mt->parse_fasta_file($aln);
	
	# iterate over genera, skip those with fewer than three species
	GENUS: for my $genus ( keys %species_for_genus ) {
		next GENUS if scalar @{ $species_for_genus{$genus} } <= 2;
	
		# aggregate sequence objects by species within this genus
		my %sequences_for_species;
		for my $species ( @{ $species_for_genus{$genus} } ) {
			if ( my @raw = $mt->get_sequences_for_taxon( $species, %fasta ) ) {
				$log->debug("$aln contained ".scalar(@raw)." sequences for species $species");
				my @seq = map { $dat->new( '-type' => 'dna', '-char' => $_ ) } @raw;
				$sequences_for_species{$species} = \@seq;
			}
		}
		
		# check if we've seen enough sequenced species
		my @species = keys %sequences_for_species;
		if ( scalar(@species) <= 2 ) {
			$log->debug("not enough species for genus $genus in $aln");
			next GENUS;
		}		
		
		# calculate the distances within the genus, take the 
		# average if species have multiple sequences
		my %dist;
		for my $i ( 0 .. ( $#species - 1 ) ) {
			for my $j ( ( $i + 1 ) .. $#species ) {
				my $sp1 = $species[$i];
				my $sp2 = $species[$j];
				$log->debug("going to compute average distance between $sp1 and $sp2");
				my @seqs1 = @{ $sequences_for_species{$sp1} };
				my @seqs2 = @{ $sequences_for_species{$sp2} };
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
		
		# pick the most distal pair, weight it by number of pairs minus one
		my ($farthest) = sort { $dist{$b} <=> $dist{$a} } keys %dist;
		$distance{$genus}->{$farthest} += scalar(keys(%dist))-1;
		$log->debug("most distal pair for genus $genus in $aln: $farthest");
	}
}

# make the final set of exemplars
my @exemplars;
for my $genus ( keys %species_for_genus ) {

	# there were 2 or fewer species in the genus, add all the species
	# from the table
	if ( not scalar keys %{ $distance{$genus} } ) {
		push @exemplars, @{ $species_for_genus{$genus} };
	}
	else {
		my %p = %{ $distance{$genus} };
		my ($sp1,$sp2) = map { split /\|/, $_ } sort { $p{$b} <=> $p{$a} } keys %p;
		push @exemplars, $sp1, $sp2;
		$log->debug(Dumper({ $genus => \%p }));
	}
}

# this includes species for which we may end up not having sequences after all
$log->info("identified ".scalar(@exemplars)." exemplars");

# make the best set of alignments:
# 1. map exemplar taxa to alignments and vice versa
my ( %taxa_for_aln, %alns_for_taxon, %nchar );
for my $aln ( @alignments ) {
	my %fasta = $mt->parse_fasta_file($aln);
	$nchar{$aln} = $mt->get_nchar(%fasta);	
	
	# iterate over all exemplars
	for my $taxon ( @exemplars ) {
		my ($seq) = $mt->get_sequences_for_taxon($taxon,%fasta);
		
		# if taxon participates in this alignment, store the mapping
		if ( $seq ) {
			$alns_for_taxon{$taxon} = [] if not $alns_for_taxon{$taxon};	
			$taxa_for_aln{$aln} = [] if not $taxa_for_aln{$aln};				
			push @{ $taxa_for_aln{$aln} }, $taxon;
			push @{ $alns_for_taxon{$taxon} }, $aln;
		}
	}
}

# 2. for each taxon, sort its alignments by decreasing taxon coverage
for my $taxon ( @exemplars ) {
	if ( my $alns = $alns_for_taxon{$taxon} ) {
		my @sorted = sort { scalar(@{$taxa_for_aln{$b}}) <=> scalar(@{$taxa_for_aln{$a}}) } @{ $alns };
		$alns_for_taxon{$taxon} = \@sorted;
	}
}

# 3. sort the taxa by increasing occurrence in alignments, so rarely sequenced species
#    are treated preferentially by including their most speciose alignments first
my @sorted_exemplars = sort { scalar(@{$alns_for_taxon{$a}}) <=> scalar(@{$alns_for_taxon{$b}}) } grep { $alns_for_taxon{$_} } @exemplars;

# starting with the least well-represented taxa...
my ( %aln, %seen );
TAXON: for my $taxon ( @sorted_exemplars ) {

	# take all its not-yet-seen alignments...
	my @alns = grep { ! $aln{$_} } @{ $alns_for_taxon{$taxon} };
	$seen{$taxon} = 0 if not defined $seen{$taxon};
	ALN: while( $seen{$taxon} < $config->BACKBONE_MIN_COVERAGE ) {
	
		# pick the most speciose alignments first
		my $aln = shift @alns;
		if ( not $aln or not -e $aln ) {
			$log->warn("no alignment available for exemplar $taxon");
			next TAXON;
		}
		$aln{$aln}++;
		
		# increment coverage count for all taxa in this alignment
		$seen{$_}++ for @{ $taxa_for_aln{$aln} };
		last ALN if not @alns;
	}
}
my @filtered_exemplars = keys %seen;
my @sorted_alignments  = keys %aln;

# produce interleaved phylip output
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
