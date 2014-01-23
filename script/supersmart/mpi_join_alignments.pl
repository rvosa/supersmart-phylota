#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util 'sum';
use Parallel::MPI::Simple;
use Bio::Phylo::Util::Logger ':levels';
use constant HEAD_NODE => 0;
use constant FASTA_FILES => 1;
use constant MATRICES => 2;

# process command line arguments
my ( $verbosity, $divergence, $overlap, $gappiness, $mergedlist, $chunkfile, $workdir ) = ( WARN, 0.20, 0.3, 8 );
GetOptions(
	'mergedlist=s' => \$mergedlist,
	'workdir=s'    => \$workdir,
	'verbose+'     => \$verbosity,
	'chunkfile=s'  => \$chunkfile,
	'divergence=f' => \$divergence,
	'overlap=f'    => \$overlap,
	'gappiness=i'  => \$gappiness,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main'
);

# here we start in parallel mode
MPI_Init();
my $rank = MPI_Comm_rank(MPI_COMM_WORLD);
	
# the following block is executed by the head node
if ( $rank == 0 ) {
		
	# dispatch jobs listed in this file
	$log->info("read alignment list file $mergedlist");
	
	# read list of file names
	open my $fh, '<', $mergedlist or die $!;	
	my @names = grep { /\S/ } <$fh>;
	chomp(@names);
	close $fh;
	        
	# MPI_Comm_size returns the total number of nodes. Because we have one
	# head node this needs to be - 1
	my $nworkers = MPI_Comm_size(MPI_COMM_WORLD) - 1;
	$log->info("we have $nworkers nodes available");
	
	# for each worker we dispatch an approximately equal chunk of clusters.
	my $nnames = int( scalar(@names) / $nworkers );
	$log->info("each worker will process about $nnames files");
	
	# this is the starting index of the subset
	my $start = 0;

	# iterate over workers to dispatch jobs
	for my $worker ( 1 .. $nworkers ) {

		# the ending index of the subset is either the starting index + chunk size or the
		# highest index in the array, whichever fits
		my $end  = ( $start + $nnames ) > $#names ? $#names : ( $start + $nnames );

		# create subset
		my @subset = @names[ $start .. $end ];
		$log->info("dispatching ".scalar(@subset)." names (index: $start .. $end) to worker $worker");
		MPI_Send(\@subset,$worker,FASTA_FILES,MPI_COMM_WORLD);

		# increment starting index
		$start += $nnames + 1;
	}
	
	# iterate over workers to receive results
	my ( @alignments, @nchars, %seqs_for_taxon );
	for my $worker ( 1 .. $nworkers ) {

		# get the result
		my $result = MPI_Recv($worker,MATRICES,MPI_COMM_WORLD);
		$log->info("received ".scalar(@{$result})." results from worker $worker");
		push @alignments, @{ $result };
	}
	
	# compute number of characters in same order as alignments
	for my $i ( 0 .. $#alignments ) {
		my $aln = $alignments[$i];
		my $nchar;
		for my $taxon ( keys %{ $aln } ) {
			$seqs_for_taxon{$taxon}++;
			if ( not defined $nchar ) {
				$nchar = length $aln->{$taxon};
			}
			elsif ( defined $nchar and $nchar != length $aln->{$taxon} ) {
				die "alignment $names[$i] is not flush!"
			}
		}
		push @nchars, $nchar;
	}	

	# read nested fractol chunks
	my ( $chunk, $ids ) = read_chunked($chunkfile);

	# this will contain two exemplars for each higher taxon
	my %exemplars;
	for my $higher_taxon ( @{ $ids } ) {
		
		# this now contains a mixture of species and higher taxa. for the
		# higher taxa we look up the exemplars.
		my @lumped = @{ $chunk->{$higher_taxon} };
	
		# this will only contain species, including the exemplars from the higher
		# taxa
		my @species;
		for my $id ( @lumped ) {
			if ( exists $exemplars{$id} ) {
				my @pair = @{ $exemplars{$id} };
				push @species, @pair;
				$log->info("expanding $id into its exemplars @pair");
			}
			else {
				push @species, $id;
			}
		}
		
		# now construct the supermatrix
		my %supermatrix;
			
		# iterate over rows
		for my $species (@species) {
			my @row;
			
			# iterate over alignments
			for my $i ( 0 .. $#alignments ) {
				
				# split seq for degap_matrix
				if ( exists $alignments[$i]->{$species} ) {
					push @row, split //, $alignments[$i]->{$species};
				}
				
				# pad with missing if species does not occur in this alignment
				else {
					push @row, '?' for 1 .. $nchars[$i];
				}
			}
			$supermatrix{$species} = \@row;
		}
		
		# remove gaps
		my ( $nchar, %degapped ) = degap_matrix(%supermatrix);
		
		# now write the supermatrix in PHYLIP format
		open my $fh, '>', "${workdir}/${higher_taxon}.phy" or die "Can't open ${workdir}/${higher_taxon}.phy : $!";
		print $fh hash_to_phylip(%degapped);
		print "${workdir}/${higher_taxon}.phy\n";
		
		# pick the two most sequence-rich exemplars
		my ( $ex1, $ex2 ) = sort { $seqs_for_taxon{$b} <=> $seqs_for_taxon{$a} } grep { $seqs_for_taxon{$_} } @species;
		$exemplars{$higher_taxon} = [ $ex1, $ex2 ];
	}	
}

# this block is executed by worker nodes
else {
	my $subset = MPI_Recv(0,FASTA_FILES,MPI_COMM_WORLD);
	$log->info("worker $rank received ".scalar(@{$subset})." files");
			
	# iterate over files
	my @result;
	for my $file ( @{ $subset } ) {
		my %matrix = process_fasta($file);		
		push @result, \%matrix if scalar keys %matrix;
	}
	
	# return result
	MPI_Send(\@result,HEAD_NODE,MATRICES,MPI_COMM_WORLD);		
}

# we need to exit cleanly by making this call on all nodes
MPI_Finalize();

sub hash_to_phylip {
	my %hash = @_;
	my $ntax    = scalar(keys(%hash));
	my ($nchar) = map { length($_) } values %hash;
	my $phylip  = "$ntax $nchar\n";
	for my $taxon ( keys %hash ) {
		$phylip .= $taxon . ( " " x ( 10 - length($taxon) ) );
		if ( ref $hash{$taxon} ) {
			$phylip .= join('', @{ $hash{$taxon} } );
		}
		else {
			$phylip .= $hash{$taxon};
		}
		$phylip .= "\n";
	}
	return $phylip;
}

sub process_fasta {
	my $fasta = shift;
	
	# parse FASTA matrix, this returns multiple sequences per taxon
	$log->info("parse FASTA file $fasta");
	my %merged = parse_matrix($fasta);		
	
	# pick the least gappy, below-threshold divergent sequence per taxon
	$log->info("pick species-level exemplars from $fasta");
	my %reduced = reduce_matrix(%merged);
	
	# remove columns with only gaps
	$log->info("remove gap-only columns from $fasta");
	my ( $nchar, %degapped ) = degap_matrix(%reduced);
	
	# count mean number of contiguous stretches. if an alignment has a
	# high number here it is because we did a profile alignment that didn't
	# work so well (perhaps actually not orthologous?) so we should omit it
	if ( scalar keys %degapped ) {
		$log->info("assess gappiness of profile alignment");
		my $gaps = 0;
		for my $key ( keys %degapped ) {
			my $row = $degapped{$key};
			my @stretches = split /-+/, $row;
			$gaps += scalar @stretches;
		}
		if ( ( $gaps / scalar keys %degapped ) > $gappiness ) {
			$log->warn("FASTA file $fasta is too gappy to include");
			%degapped = ();
		}
	}

	return %degapped;
}

# removes columns with only '-' and/or '?' characters
sub degap_matrix {
	my %matrix = @_;
	my %degapped;
	$log->info("remove columns with only gaps");
	
	my $i = 0;
	my $nchar = 0;
	my $ntax = scalar(keys(%matrix));
	DEGAP: while(1) {
		last DEGAP if $ntax == 0;
		
		# detect gaps at column $i
		my $columns_with_data = 0;
		for my $id ( keys %matrix ) {
			
			# this is how we break out the infinite loop
			last DEGAP if not exists $matrix{$id}->[$i];
			$degapped{$id} = '' if not exists $degapped{$id};
			my $token = $matrix{$id}->[$i];
			
			# the flag is switched when at least one token is not - or ?
			if ( $token ne '-' and $token ne '?' ) {
				$columns_with_data++;
			}
		}
		
		# skip or grow degapped matrix
		if ( ( $columns_with_data / $ntax ) < $overlap ) {
			$log->debug("not enough data at column $i");
		}
		else {
			for my $id ( keys %matrix ) {
				$degapped{$id} .= $matrix{$id}->[$i];
				$nchar = length($degapped{$id});
			}
		}
		$i++;
	}
	
	# remove rows with only ???
	my @rows = keys %degapped;
	for my $row ( @rows ) {
		delete $degapped{$row} if $degapped{$row} =~ /^\?+$/;
	}
	
	return $nchar, %degapped;
}

# reads a fractol chunk file
sub read_chunked {
	my $infile = shift;
	my ( %chunk, @id );
	$log->info("read chunk file $infile");
	
	# open file handle
	open my $fh, '<', $infile or die $!;
	
	# iterate over lines
	LINE: while(<$fh>) {
		chomp; # strip line ending
		next LINE unless /\t/;
		
		# split the line on tab character into key and value
		my ( $k, $v ) = split /\t/, $_;
		push @id, $k;
		
		# split the value on commas, assign to hash
		$chunk{$k} = [ split /,/, $v ];
	}
	close $fh;
	return \%chunk, \@id;
}

# reads a FASTA alignment, returns hash keyed on taxon ID
sub parse_matrix {
	my $file = shift;
	my ( %matrix, $current );
	
	# open file handle 
	open my $fh, '<', $file or die $!;
	$log->info("read sequences from FASTA file $file");
	
	# read over the file handle
	while(<$fh>) {
		chomp; # strip line ending
		
		# this matches the FASTA definition line, we capture the taxon ID
		if ( />.*taxon\|(\d+)/ ) {
			$current = $1;
			$log->debug("found taxon ID $current");
			
			# already seen this ID, now seeing it for the 2nd (or more) time.
			if ( $matrix{$current} ) {
				$log->debug("already seen $current, starting new empty string");
				push @{ $matrix{$current} }, '';
			}
			
			# this is the first time we see the taxon ID
			else {
				$log->debug("not yet seen $current, initializing array with an empty string");
				$matrix{$current} = [ '' ];
			}
		}
		else {			
			s/\s//g;
			$matrix{$current}->[-1] .= $_; # append to last seq
		}
	}
	return %matrix;
}

# compute uncorrected distance between two sequences
sub compute_dist {
	my ( $seq1, $seq2 ) = @_;
	
	# the focal sequence, one site per element
	my @s1 = split //, $seq1;
	
	# comparison sequence, one site per element
	my @s2 = split //, $seq2;
	
	# pairwise distance
	my $dist = 0;
			
	# compare all sites
	for my $i ( 0 .. $#s1 ) {
		$dist++ if $s1[$i] ne $s2[$i] and $s1[$i] ne '-' and $s2[$i] ne '-';
	}
			
	# divide by seq length
	$dist /= scalar(@s1);	
	return $dist;
}

# picks the least gappy, below-threshold divergent sequence per taxon
sub reduce_matrix {
	my %matrix = @_;
	
	# this will be a hash where the key is the raw seq, the value are the
	# pairwise distances with all other sequences	
	my %distance = map { $_ => [] } map { @{ $_ } } values %matrix;
	my @seqs = keys %distance;
	$log->info("compute all pairwise sequence distances");
	
	# because they are pairs (i.e. different sequences), the outer loop has to
	# end at $#seqs - 1 and the inner loop has to start at $i + 1
	for my $i ( 0 .. ( $#seqs - 1 ) ) {
		$log->debug("progress: ".int($i/$#seqs*100)."%");
		
		# here we iterate over all pairs
		for my $j ( ( $i + 1 ) .. $#seqs ) {
			
			# add pairwise distance to tallies
			my $dist = compute_dist( $seqs[$i], $seqs[$j] );
			push @{ $distance{$seqs[$i]} }, $dist;
			push @{ $distance{$seqs[$j]} }, $dist;
			$log->debug("distance $i <-> $j is $dist");			
		}		
	}
	
	# now compute the averages
	for my $seq ( keys %distance ) {
		my @distances = @{ $distance{$seq} };
		my $mean = sum(@distances) / scalar(@distances);
		$distance{$seq} = $mean;
		if ( $mean > $divergence ) {
			$log->warn("detected highly divergent sequence: $mean");
		}
	}
	
	# now we're going to reduce the sequence set to use only the longest
	# sequence in the set for that taxon ID
	$log->info("select longest sequence from among within-species sequences");
	for my $taxon_id ( keys %matrix ) {
		
		# filter out too-divergent sequences
		my @sequences = grep { $distance{$_} < $divergence } @{ $matrix{$taxon_id} };
		
		# mmm...
		if ( not @sequences ) {
			$log->warn("all sequences for $taxon_id were too divergent");
			delete $matrix{$taxon_id};
		}
		else {
			# we sort in ascending order from fewer to more gaps
			my %gapcount;
			for my $s ( @sequences ) {
				my @regions = split /-+/, $s;
				$gapcount{$s} = scalar @regions;
			}
			my @sorted = sort { $gapcount{$a} <=> $gapcount{$b} } @sequences;
			$matrix{$taxon_id} = [ split //, $sorted[0] ];
		}
	}
	return %matrix;
}
