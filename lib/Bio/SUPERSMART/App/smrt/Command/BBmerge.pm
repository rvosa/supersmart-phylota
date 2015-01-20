package Bio::SUPERSMART::App::smrt::Command::BBmerge;

use strict;
use warnings;

use List::Util qw(min max sum);
use Bio::Phylo::Matrices::Datum;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::TreeService;

use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: creates supermatrix for genus-level backbone tree

=head1 NAME

BBmerge.pm - creates supermatrix for genus-level backbone tree

=head1 SYNOPSIS

smrt bbmerge [-h ] [-v ] [-w <dir>] -a <file> -t <file> [-o <file>] 

=head1 DESCRIPTION

Given an input file that lists alignment file locations (one on each line), traverses each 
genus in each alignment and picks the most divergent two species to represent 
their genus. The rationale is that these species will (likely) cross the root of their 
genus, so that the below genus-level tree can then be scaled to the same depth of that 
split and be grafted onto the tree without (too many) negative branch lengths.

The way in which the overall two divergent species within the genus are selected is as 
follows:

=over

=item * for each alignment, within each genus, make all pairwise comparisons and sort
the pairs by decreasing sequence divergence.

=item * pick the most distal pair and weight it in proportion to the number of pairs, 
within that genus for that alignment, minus one. This means that singleton pairs are 
discarded, and those from bigger samples are assumed to more accurately indicate which 
taxa actually cross the root.

=item * after having processed all alignments, pick the species pair that has the highest
score. 

=back

Subsequently, the optimal combination of markers needs to be selected to best cover the
exemplars. It is not optimal to just concatenate all alignments that cover any of the 
taxa - this can result in monstrous, sparse, supermatrices. Instead we give the user the
possibility of assembling a set of alignments such that all exemplar species are covered
by at least some minimal value, (though relatively frequently studied species would exceed
this). This is done as follows:

=over

=item * for each exemplar species, collect all alignments that include it and sort this
collection in decreasing exemplar taxon coverage (i.e. the first alignment has the most
exemplar species in it, the last alignment the fewest).

=item * sort the exemplar species by increasing overall participation in the alignments
(i.e. the first exemplar has been sequenced the fewest times, the last one the most).

=item * iterate over the sorted list of exemplars, and for each exemplar add their 
not-yet-seen, sorted alignments to the stack, one by one. After adding each alignment,
update the coverage counts for all exemplar species that participate in that alignment.
End the iterations when all exemplars have crossed their threshold or have no more 
alignments available.

=back

B<Point of consideration>: the node depths on the exemplar tree will be underestimates
relative to the genus-level tree (due to the node density effect), so it might be better
to give the node depths from the exemplar tree as priors and then see where the posteriors
come out. Disadvantage is that this is likely to lead to negative branch lengths in some
cases.

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "supermatrix.phy";
	my $outformat_default = "phylip";
	my $markerstable_default = "markers-backbone.tsv";
	return (
		["alnfile|a=s", "list of file locations of merged alignments  as produced by 'smrt orthologize'", { arg => "file", mandatory => 1}],	
		["taxafile|t=s", "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'", { arg => "file", mandatory => 1}],
		["outfile|o=s", "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],	
		["format|f=s", "format of supermatrix, defaults to '$outformat_default'; possible formats: bl2seq, clustalw, emboss, fasta, maf, mase, mega, meme, msf, nexus, pfam, phylip, prodom, psi, selex, stockholm", {default => $outformat_default}],	
		["markersfile|o=s", "name for summary table with included accessions, defaults to $markerstable_default", { default=> $markerstable_default, arg => "file"}],
	
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	#  If alignment or taxa file is absent or empty, abort  
	my @files = ( $opt->alnfile, $opt->taxafile );
	foreach my $file ( @files ){
		$self->usage_error("need alnfile and taxafile arguments") if not $file;
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
}


sub run {
	my ($self, $opt, $args) = @_;		
		
	# collect command-line arguments
	my $alnfile = $opt->alnfile;
	my $taxafile = $opt->taxafile;
	my $outfile= $self->outfile;	
	my $format = $opt->format;
	my $markersfile  = $opt->markersfile;

	# instantiate helper objects
	my $mts     = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $ts		= Bio::Phylo::PhyLoTA::Service::TreeService->new;	
	my $dat    	= 'Bio::Phylo::Matrices::Datum';
	my $config 	= Bio::Phylo::PhyLoTA::Config->new;
	my $mt     	= Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $log = $self->logger;

	# read list of alignments
	my @alignments;

	$log->info("going to read alignment list $alnfile");
	open my $fh, '<', $alnfile or die $!;
	while(<$fh>) {
		chomp;
		push @alignments, $_ if /\S/ && -e $_;
	}
	
	# read taxa table
	$log->info("going to read taxa table $taxafile");
	my @records = $mt->parse_taxa_file($taxafile);
	
	# extract the distinct species for each genus
	my %species_for_genus;
	
	my @g = $mt->get_distinct_taxa( 'genus' => @records );
	
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
			my $genus_name = $ts->find_node($genus)->taxon_name;
			if ( scalar @{ $species_for_genus{$genus} } <= 2 ) {
				$log->debug( "Genus $genus ($genus_name) has less than three species. Skipping." );
				next GENUS; 
			}
		
			# aggregate sequence objects by species within this genus
			my %sequences_for_species;
			for my $species ( @{ $species_for_genus{$genus} } ) {
				my $species_name = $ts->find_node($species)->taxon_name;
				if ( my @raw = $mt->get_sequences_for_taxon( $species, %fasta ) ) {
					$log->debug("$aln contained ".scalar(@raw)." sequences for species $species ($species_name)");
					my @seq = map { $dat->new( '-type' => 'dna', '-char' => $_ ) } @raw;
					$sequences_for_species{$species} = \@seq;
				}
			}
			
			# check if we've seen enough sequenced species
			my @species = keys %sequences_for_species;
			if ( scalar(@species) <= 2 ) {
				$log->debug("not enough sequenced species for genus $genus ($genus_name) in $aln");
				next GENUS;
			}		
			
			# calculate the distances within the genus, take the 
			# average if species have multiple sequences
			my %dist;
			for my $i ( 0 .. ( $#species - 1 ) ) {
				for my $j ( ( $i + 1 ) .. $#species ) {
					my $sp1 = $species[$i];
					my $sp2 = $species[$j];
					#$log->debug("going to compute average distance between $sp1 and $sp2");
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
			##$log->debug("most distal pair for genus $genus in $aln: $farthest");
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
			my ($sp1,$sp2) = map { split (/\|/, $_) } sort { $p{$b} <=> $p{$a} } keys %p;
			push @exemplars, $sp1, $sp2;
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
	    my $taxon_name = $ts->find_node($taxon)->taxon_name;
		$log->info("Checking alignment coverage for taxon with name $taxon_name");
		# take all its not-yet-seen alignments...
		my @alns = grep { ! $aln{$_} } @{ $alns_for_taxon{$taxon} };
		$seen{$taxon} = 0 if not defined $seen{$taxon};
		ALN: while( $seen{$taxon} < $config->BACKBONE_MIN_COVERAGE ) {
			# pick the most speciose alignments first
			my $aln = shift @alns;
	                if ( not $aln or not -e $aln ) {
						$log->warn("no alignment available for exemplar $taxon ($taxon_name)");
						next TAXON;
			}
			$aln{$aln}++;
			
			# increment coverage count for all taxa in this alignment
			$seen{$_}++ for @{ $taxa_for_aln{$aln} };
			last ALN if not @alns;
		}
	}
	
	# filter exemplars: only take the ones with sufficient coverage
	my @filtered_exemplars = grep { $seen{$_} >= $config->BACKBONE_MIN_COVERAGE } keys %seen;
	my @sorted_alignments  = keys %aln;
	
	# filter the alignments: only include alignments in supermatrix which have
	#  at least one species from the filtered exemplars ==> no empty rows in supermatrix
	my @filtered_alignments;
	for my $aln ( @sorted_alignments ) {
	        my $count = 0;
	        # increse counter if an alignment is present for a taxon
	        for my $taxon ( @filtered_exemplars ) {
	                my %fasta = $mt->parse_fasta_file($aln);
	                my ( $seq ) = $mt->get_sequences_for_taxon($taxon,%fasta);
					if ( $seq ) {
			                        $count++;
			                        last;
					}
	        }
	        if ( $count > 0 ) {
	                $log->info("including $aln in supermatrix");
	                push @filtered_alignments, $aln; 
	                #my @markers = get_seed_gi_description($aln);
	                #$log->info("marker : ".$markers[0])
	        } 
	        else {
	                $log->info("alignment $aln not included in supermatrix");
	        }
	}
	
	$log->info("using ".scalar(@filtered_alignments)." alignments for supermatrix");
	$log->info("number of filtered exemplars : " . scalar(@filtered_exemplars));
	
	# filter supermatrix for gap-only columns and write to file
	$self->_write_supermatrix( \@filtered_alignments, \@filtered_exemplars, $outfile, $format, $markersfile );
	
	$log->info("DONE, results written to $outfile");
	
	return 1;
	
}

sub _write_supermatrix {
	my ($self, $alnfiles, $exemplars, $filename, $format, $markersfile) = @_;
	my $log = $self->logger;

	# make has with concatenated sequences per exemplar
	my %allseqs = map{ $_ => "" } @{$exemplars}; 
	
	# also store which markers have been chosen for each exemplar
	#my %allmarkers;
	my @marker_table;
	
	my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	
	# Retrieve sequences for all exemplar species
	for my $alnfile ( @$alnfiles ) {
		my %fasta = $mt->parse_fasta_file($alnfile);
		my $nchar = length $fasta{(keys(%fasta))[0]};
		my %row;
		for my $taxid ( @$exemplars ) {			
			my ( $seq ) = $mt->get_sequences_for_taxon($taxid, %fasta);
			
			my ($header) = grep {/$taxid/} keys(%fasta);
			$row{$taxid} = $header if $seq;
			
			if ( ! $seq ) {	
	        	$seq = '?' x $nchar;
			}
			$allseqs{$taxid} .= $seq;
		}
		push @marker_table, \%row;
	}


	$mts->write_marker_summary( $markersfile, \@marker_table, $exemplars );


	# Delete columns that only consist of gaps
	my $nchar = length $allseqs{(keys(%allseqs))[0]};
	my $vals = values %allseqs;
	my %ind;
	for my $v (values %allseqs) {
	  # check if column only consists of gap characters
	  $ind{ pos($v) -1 }++ while $v =~ /[-Nn\?]/g;
	}
	my @to_remove = sort {$b <=> $a} grep { $ind{$_} == $vals } keys %ind;
	
	# remove columns at specified indices
	for my $v (values %allseqs) {
	  substr($v, $_, 1, "") for @to_remove;
	}

	my $removed = $nchar - length $allseqs{(keys(%allseqs))[0]};
	$log->info("Removed $removed gap-only columns from supermatrix");

	# Write supermatrix to file
	my $aln = Bio::SimpleAlign->new();
	map {$aln->add_seq(Bio::LocatableSeq->new(-id=>$_, seq=>$allseqs{$_}, start=>1)) } keys %allseqs;
	
	my $stream = Bio::AlignIO->new(-format	 => $format,
                                   -file     => ">$filename",
                                   -idlength => 10 );
	
	$stream->write_aln($aln);	
}


1;