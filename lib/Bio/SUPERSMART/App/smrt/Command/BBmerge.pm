package Bio::SUPERSMART::App::smrt::Command::BBmerge;

use strict;
use warnings;

use Bio::Phylo::Matrices::Datum;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::SUPERSMART::App::SubCommand;

use List::MoreUtils qw(uniq);
use List::Util qw(max);

use base 'Bio::SUPERSMART::App::SubCommand';
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
        ["include_taxa|i=s", "one or multiple names of taxa present in <taxafile> (e.g. species or genus names, separated by commata) whose representative species will be included in the output dataset, regardless of marker coverage and sequence divergence", {}],
		["markersfile|m=s", "name for summary table with included accessions, defaults to $markerstable_default", { default=> $markerstable_default, arg => "file"}],
	
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
	my $include_taxa = $opt->include_taxa;

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
		
	# valid ranks for taxa in supermatrix                                                                                             
	my @ranks = ("species", "subspecies", "varietas", "forma");
	
	# remember species (or lower) for taxa that have to be included in supermatrix, if given as argument
	my %user_taxa;
	if ( $include_taxa ) {
	    my @taxa = split(',', $include_taxa);
	    my @ids = map { my @nodes = $ts->search_node( {taxon_name=>$_} )->all; $nodes[0]->ti } @taxa;
	    my @species = $mt->query_taxa_table(\@ids, \@ranks, @records);
	    %user_taxa = map { $_ => 1} @species;
	}

 	# extract the distinct species for each genus
	my %species_for_genus;
	
	my @g = $mt->get_distinct_taxa( 'genus' => @records );
	 	
	for my $genus ( $mt->get_distinct_taxa( 'genus' => @records ) ) {
	    # extract the distinct species (and lower) for the focal genus                                                                                                                                                                                                  
	    my @species = $mt->query_taxa_table( $genus, \@ranks, @records );
	    $species_for_genus{$genus} = \@species;
	}
    
	$log->info("read ".scalar(keys(%species_for_genus))." genera");
	
        	
	my @all_species = map {@$_} values(%species_for_genus);	

	# build an adjacency matrix connecting species if they share markers. 
	my %adjacency_matrix = map {$_=> {map{$_=>0} @all_species}} @all_species;

	for my $aln ( @alignments ) {
	    my %fasta = $mt->parse_fasta_file($aln);
	    my @species = uniq map {$1 if $_=~m/taxon\|([0-9]+)/} keys(%fasta);
	    for my $s1 (@species) {
		for my $s2 (@species) {
		    $adjacency_matrix{$s1}{$s2}++;
		}
	    }	    
	}
	# do not keep the diagonal
	for my $sp (keys %adjacency_matrix) {
	    delete $adjacency_matrix{$sp}->{$sp};
	}
	
	# now pick the exemplars: rather than just taking the two most distal species
	#  within a genus, we narrow down the list of possible exemplars by two criteria 
	#  first: species must share markers with species in other genera, and species
	#  must share markers with species in it's own genus. After that, we take the most
	#  distal pair of the remaining species.	
	my @exemplars;
	my %distance = map { $_ => {} } keys %species_for_genus;	

     GENUS: for my $genus  ( keys %species_for_genus ) {
	  my $genus_name = $ts->find_node($genus)->taxon_name;	  
	  my %species = map {$_=>1} @{$species_for_genus{$genus}};

	  if ( scalar(keys(%species)) <= 2 ) {
	      $log->debug( "Genus $genus ($genus_name) has less than three species. Skipping." );
	      next GENUS; 
	  }

	  my @candidates;
	  for my $s (keys(%species)) {
	      my %adj = %{$adjacency_matrix{$s}};
	      # species is a candidate to be an exemplar if it is connected to species 
	      #  inside -and- outside it's own genus
	      my @connected_genus = grep { exists $species{$_} and $adj{$_} > 0 } keys(%adj);
	      my @connected_other = grep { ! exists $species{$_} and $adj{$_} > 0 } keys(%adj);
	      push @candidates, $s if scalar(@connected_genus) and scalar (@connected_other);	      
	  }
	  
	  # now calculate the distance between the candidate species for the genus
	  my %distance;
	ALN: for my $aln ( @alignments ) {
	    $log->debug("checking for species in $aln");
	    my %fasta = $mt->parse_fasta_file($aln);
		    
	    # aggregate sequence objects by species within this genus	      
	    my %sequences_for_species;
	    for my $candidate ( @candidates ) {
		my $name = $ts->find_node($candidate)->taxon_name;
		if ( my @raw = $mt->get_sequences_for_taxon( $candidate, %fasta ) ) {
		    $log->debug("$aln contained ".scalar(@raw)." sequences for species $candidate ($name)");
		    my @seq = map { $dat->new( '-type' => 'dna', '-char' => $_ ) } @raw;
		      $sequences_for_species{$candidate} = \@seq;
		}
	    }
	    
	    # check if we've seen enough sequenced species	     
	    if ( scalar(keys(%sequences_for_species)) < 3 ) {
		$log->debug("not enough sequenced species for genus $genus ($genus_name) in $aln");
		next ALN;
	    }
	    # calculate the distances within the genus, take the 
	    # average if species have multiple sequences
	    my %dist;
	    for my $i ( 0 .. ( $#candidates - 1 ) ) {
		for my $j ( ( $i + 1 ) .. $#candidates ) {
		    my $sp1 = $candidates[$i];
		    my $sp2 = $candidates[$j];
		    if ( $sequences_for_species{$sp1} and $sequences_for_species{$sp2} ) {
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
	    }
	      # pick the most distal pair, weight it by number of pairs minus one
	    my ($farthest) = sort { $dist{$b} <=> $dist{$a} } keys %dist;
	    $distance{$farthest} += scalar(keys(%dist))-1;
	}
	  if ( not scalar keys %distance ) {
	      push @exemplars, @candidates;
	      my @names = map {$ts->find_node($_)->taxon_name} @candidates ;
	      $log->debug('adding species ' . join(', ', @names) . ' to exemplars');		
	  }
	  else {
	      my ($sp1,$sp2) = map { split (/\|/, $_) } sort { $distance{$b} <=> $distance{$a} } keys %distance;
	      push @exemplars, $sp1, $sp2;
	      $log->debug('adding species ' . $ts->find_node($sp1)->taxon_name . ', ' . $ts->find_node($sp2)->taxon_name . ' to exemplars');
	  }	  
     }
	use Data::Dumper;
	print Dumper(\@exemplars);
	die;
	
	  
	# make the final set of exemplars
	#my @exemplars;
	for my $genus ( keys %species_for_genus ) {
	    
	    # there were 2 or fewer species in the genus, add all the species
	    # from the table
	    if ( not scalar keys %{ $distance{$genus} } ) {
		push @exemplars, @{ $species_for_genus{$genus} };
		my @names = map {$ts->find_node($_)->taxon_name} @{ $species_for_genus{$genus} };
		$log->debug('adding species ' . join(', ', @names) . ' to exemplars');		
	    }
	    else {
		my %p = %{ $distance{$genus} };
		my ($sp1,$sp2) = map { split (/\|/, $_) } sort { $p{$b} <=> $p{$a} } keys %p;
		push @exemplars, $sp1, $sp2;
		$log->debug('adding species ' . $ts->find_node($sp1)->taxon_name . ', ' . $ts->find_node($sp2)->taxon_name . ' to exemplars');
	    }
	}
	
	# This includes species for which we may end up not having sequences after all
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
	my @filtered_exemplars = grep { $seen{$_} >= $config->BACKBONE_MIN_COVERAGE or exists $user_taxa{$_} } keys %seen;
	
	my @sorted_alignments  = sort keys %aln;
	
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
	        } 
	        else {
	                $log->info("alignment $aln not included in supermatrix");
	        }
	}
	
	$log->info("using ".scalar(@filtered_alignments)." alignments for supermatrix");
	$log->info("number of filtered exemplars : " . scalar(@filtered_exemplars));
	
	# filter supermatrix for gap-only columns and write to file
	$self->_write_supermatrix( \@filtered_alignments, \@filtered_exemplars, \%user_taxa, $outfile, $format, $markersfile );
	
	$log->info("DONE, results written to $outfile");
	
	return 1;
	
}

sub _write_supermatrix {
	my ($self, $alnfiles, $exemplars, $user_taxa, $filename, $format, $markersfile) = @_;
	
	my $log = $self->logger;
    my %user_taxa = %{$user_taxa};

	# make has with concatenated sequences per exemplar
	my %allseqs = map{ $_ => "" } @{$exemplars}; 
	
	# also store which markers have been chosen for each exemplar
	my @marker_table;
	
	my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	
	# Retrieve sequences for all exemplar species and populate marker table
	for my $alnfile ( @$alnfiles ) {
		my %fasta = $mt->parse_fasta_file($alnfile);
		my $nchar = length $fasta{(keys(%fasta))[0]};
		my %row;
		for my $taxid ( @$exemplars ) {			
			my @seqs = $mt->get_sequences_for_taxon($taxid, %fasta);
			$log->warn("Found more than one sequence for taxon $taxid in alignment file $alnfile. Using first sequence.") if scalar(@seqs) > 1; 
			my $seq;
			if ( $seq = $seqs[0] ) {			
				my ($header) = grep {/$taxid/} sort keys(%fasta);
				$row{$taxid} = $header;
			}
			else {	
	        	$seq = '?' x $nchar;
			}
			$allseqs{$taxid} .= $seq;
		}
		push @marker_table, \%row;
	}

	$log->info("Filtering exemplars for subsets that do not have marker overlap");
	my @pruned_exemplars = $self->_prune_exemplars (\@marker_table);
		
	
	my @excluded = grep {my $ex = $_; ! grep($_ == $ex, @pruned_exemplars) } @$exemplars;
	$log->warn("Found "  . scalar(@excluded) . " exemplars that do not share markers with other exemplars") if scalar (@excluded);
	
    # if user specified exemplars were removed, add them back!                        
    my @excluded_user_taxa = grep {exists $user_taxa{$_}} @excluded;                                                                                        

    push @pruned_exemplars, @excluded_user_taxa;                                                                                   

    foreach my $t (@excluded_user_taxa) {                                       
           $log->warn("Taxon $t should be excluded due to insufficient marker coverage but is given by in include_taxa argument");                                              
    }  
	
	foreach my $exc (@excluded){
		my $name = $mts->find_node($exc)->get_name;
		$log->warn("Excluding exemplar with name $name and id $exc");		
	}
		
	# remove rejected exemplars from marker table
	my %valid_exemplars = map {$_=>1} @pruned_exemplars;
	foreach my $m (@marker_table){		
		foreach my $k ( keys %{$m} ) {			
			if ( ! exists $valid_exemplars{$k} ){
				delete $m->{$k};
			}
		}		
	}
	
	# write table listing all marker accessions for taxa
	$mts->write_marker_summary( $markersfile, \@marker_table, \@pruned_exemplars );
	
	# only keep the sequences that belong to the pruned exemplars
	%allseqs = map { $_=> $allseqs{$_} } grep {exists $allseqs{$_} } @pruned_exemplars;
		
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


# It often happens that there are subsets of exemplars that do not share markers with other sets of exemplars.
#  When inferring a tree with such a disconnected supermatrix, disconnected sets show very long branch lengths.
#  This subroutine returns the largest subset of exemplars that have common markers. This is done by doing a 
#  breadth-first search on the adjaceny matrix of exemplars connected by their respective markers.
sub _prune_exemplars {
	my ($self, $tab ) = @_;
	my @marker_table = @$tab;
	
	# create adjacency matrix
	my %adj;
	foreach my $column ( @marker_table ) {
		my @species = sort keys (%$column); 
		for my $i ( 0..$#species ) {
			
			for my $j ( $i+1..$#species ) {				
				my %specs_i = map { $_ => 1 } @{$adj{$species[$i]}};
				my %specs_j = map { $_ => 1 } @{$adj{$species[$j]}};
				
				if ( ! exists $specs_i{$species[$j]}) {
					push @{$adj{$species[$i]}}, $species[$j];			
				}
				if ( ! exists $specs_j{$species[$i]}) {								
					push @{$adj{$species[$j]}}, $species[$i];
				}							
			}			
		}						
	}
	
	# do a BFS to get the unconnected subgraphs from the adjacency matrix	
	my @sets;	
	my @current_set;
	
	my @queue = ( sort keys (%adj))[0];
		
	while ( @queue ) {
		my $current_node = shift (@queue);						
		push @current_set, $current_node;				
		if (exists $adj{$current_node}){						
			my @neighbors = @{$adj{$current_node}};
			if ( @neighbors ) {
				foreach my $node ( @neighbors ) {
					if ($node != $current_node){
						push @queue, $node;			
					}		
				}		
			}
			delete $adj{$current_node};
		}				
		@current_set = uniq (@current_set);
		@queue = uniq (@queue);

		if (scalar (@queue) == 0){
			my @cs = @current_set;
			push @sets, \@cs;			
			@current_set = ();			
			if ( keys %adj ){
				my ($key) = sort keys %adj; 
				push @queue, $key;
			}
		}		
	}	
	
	my $max_exemplars = max( map {scalar (@$_)} @sets);
	foreach my $s (@sets) {
		if ( scalar @$s == $max_exemplars ) {
			return @$s;
		}
		
	}	
}


1;
