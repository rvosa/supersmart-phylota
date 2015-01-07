package Bio::SUPERSMART::App::smrt::Command::Orthologize;

use strict;
use warnings;

use File::Copy qw(copy);
	
use Bio::SearchIO;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;

use Bio::Phylo::PhyLoTA::Service::ParallelService 'pthreads';

use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: creates orthologous clusters of aligned sequences

=head1 NAME

Align.pm - assesses orthology in different sequence alignments and merges them into orthologous clusters

=head1 SYNOPSYS


=head1 DESCRIPTION

Given a list of aligned candidate clusters, assigns orthology among the clusters by performing reciprocal 
blast searches on the seed sequences around which the clusters were assembled. Produces a list of 
re-aligned superclusters. 

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "merged.txt";
	return (
		["infile|i=s", "list of file locations of multiple sequence alignments  as produced by 'smrt align'", { arg => "file", mandatory => 1}],
		["outfile|o=s", "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => "filename"}],	
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	# We only have to check the 'infile' argument. 
	#  If the infile is absent or empty, abort  
	my $file = $opt->infile;
	$self->usage_error("no infile argument given") if not $file;
	$self->usage_error("file $file does not exist") unless (-e $file);
	$self->usage_error("file $file is empty") unless (-s $file);			
}


sub run {
	my ($self, $opt, $args) = @_;
	
	# collect command-line arguments
	my $infile = $opt->infile;
	my $outfile = $opt->outfile;
	(my $workdir = $opt->workdir) =~ s/\/$//g;
	
	# create empty output file 
	open my $outfh, '>', $workdir . '/' . $outfile or die $!;
	close $outfh;
	
	# instantiate helper objects
	my $service = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
	my $mts     = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $config  = Bio::Phylo::PhyLoTA::Config->new;
	my $log     = $self->logger;
	
	my $dbname = File::Spec->catfile( $workdir, 'seeds.fa' );
	$log->info("creating database file $dbname");
	
	# read list of alignments
	$log->info("going to read file list $infile");
	open my $fh, '<', $infile or die $!;
	my @alignments = <$fh>;
	chomp @alignments; 
	close $fh;
		
	# write seed sequences, filter out seqs that have the same seed
	my %seen;
	$log->info("going to write sequences to $dbname");
	open my $sfh, '>', $dbname or die $!;
	for my $aln ( @alignments ) {
		if ( $aln =~ /(\d+)-.+\.fa/ ) {
			my $gi = $1;
			if ( ! exists($seen{$gi}) )  {
				$seen{$gi}++;
				my $seq = $service->find_seq($gi);
				print $sfh '>', $gi, "\n", $seq->seq, "\n";
				print $sfh "alignment\n";
			}
		}	
	}
	$log->info("read ".scalar(@alignments)." alignments");

	# make blast db
	$log->info("going to make BLAST db for $dbname");
	my @cmd = (
		$config->MAKEBLASTDB_BIN,
		'-in'     => $dbname,
		'-dbtype' => 'nucl',
		'2>'      => '/dev/null',
		'>'       => '/dev/null',		
	); 	
	$log->info( "going to run command ".join(" ", @cmd) );
	my $result = `@cmd`;
	
	# run blast search
	$log->info("going to run all vs all BLAST search on $dbname");
	@cmd = (
		$config->BLASTN_BIN,
		'-query'  => $dbname,
		'-db'     => $dbname,
		'2>'      => '/dev/null',	
	);
	$log->info( "going to run command ".join(" ", @cmd) );
	$result = `@cmd`;
	
	die "BLAST result empty" if not length($result) > 0;
		
	$log->info("going to process BLAST results");	
	
	open my $out, '<', \$result;

	my $report = Bio::SearchIO->new( '-format' => 'blast', '-fh' => $out );
	
	my @blastresults;
	while ( my $result = $report->next_result ) {
		push @blastresults, $result;
	}	
	$log->info("number of blast results : ".scalar(@blastresults));
			
	# process results, this is all-vs-all so many results with many hits. we will discard
	# all hits where the overlapping region is smaller than 0.51 (default) times of the length 
	# of both query and hit
	my $overlap = $config->MERGE_OVERLAP;
	my @res = pmap {
		my ($result) = @_;
		my $query = $result->query_name;
		$log->info("querying for seed gi $query");
		my $q_l = length($service->find_seq($query)->seq);
		my %blhits;
		$blhits{$query} = [];
		
		# iterate over hits for focal query
		while ( my $hit = $result->next_hit ) {
			if ( my $name = $hit->name ) {
				my $h_l = length($service->find_seq($name)->seq);
				
				# add up the lengths of the overlapping regions
				my $q_hsp_l = 0;
				my $h_hsp_l = 0;			
				while( my $hsp = $hit->next_hsp ) {
					$q_hsp_l += $hsp->length('query');
					$h_hsp_l += $hsp->length('hit');
				}
				
				# check if the overlapping regions are long enough
				if ( ($q_hsp_l/$q_l) > $overlap and ($h_hsp_l/$h_l) > $overlap ) {
					push @{ $blhits{$query} }, $name;
					$log->debug("\thit: $name");
				}
				else {
					$log->debug("discarding hit $name for query $query");
				}
			}
		}
		
		# report progress
		my $template = 'found %i hits for query %s (%i nt)';
		my @args = ( scalar(@{$blhits{$query}}), $query, $q_l );
		my $message  = sprintf $template, @args;
		$log->info($message);	
		return \%blhits;
	} @blastresults;

	my %hits;
	map {@hits{keys %{$_} } = values %{$_} if $_} @res;
	
	# make single linkage clusters
	my $sets = [];
	for my $gi ( keys %hits ) {
		if ( $hits{$gi} ) {
			$log->info("going to cluster from seed $gi");
			_cluster( $sets, delete $hits{$gi}, \%hits);
		}
	}

	# now remove duplicates
	@$sets = values %{ { map { join('|', sort { $a <=> $b } @{$_}) => $_ } @$sets } };
	
	# assign ids to clusters
	my $clustid = 1;
	my @clusters = map { {'id'=>$clustid++, 'seq'=>$_} } @$sets;
	
	# align
	my $total = scalar @clusters;

	my @cres = pmap {
		my ($clref) = @_; 
	    my %cluster = %{$clref};
		my $clusterid = $cluster{'id'};
		my @seqids = @{$cluster{'seq'}};
					
		# turn GIs into file names 
		my @files = map { glob ( "$workdir/" .  $_. "*.fa" ) } @seqids;
	
		# the name of the file that contains the merger of this cluster
		my $merged = File::Spec->catfile( $workdir, "cluster${clusterid}.fa" );
		
		# check for singletons
		if ( scalar(@files) == 1 ) {
			copy $files[0], $merged;
			$log->info("singleton cluster $clusterid: @files");
			open my $outfh, '>>', $workdir . '/' . $outfile or die $!;
			print $outfh $merged, "\n";
			close $outfh;
			return (0);
		} 
		# (the effectiveness of this procedure may depend on the input order). I was trying
		# to do this using List::Util::reduce, but it segfaults.
		$log->info("going to reduce @files");		
				
		if ( $files[0] ) {
			copy $files[0], $merged;
		}
				
		for my $i ( 1 .. ( $#files ) ) {
			# as long as merges are unsuccessful, we will continue to attempt
			# profile align the subsequence files against the first input file. TODO:
			# sort the input files by size so that if all mergers fail, we end up with
			# the largest input file

			my $file2 = $files[$i];
						
			# do the profile alignment
			$log->info("attempting to merge $merged and $file2 (# " .  $i  . " of " . $#files . ")");
			
			if ( my $result = $service->profile_align_files($merged,$file2) ){
				# evaluate how this went
				if ( $mts->calc_mean_distance($result) < $config->BACKBONE_MAX_DISTANCE ) {
					my %fasta  = $mts->parse_fasta_string($result);				
					%fasta = $mts->dedup(%fasta);
					open my $fh, '>', $merged or die $!;
					for my $defline ( keys %fasta ) {
						print $fh '>', $defline, "\n", $fasta{$defline}, "\n";
					}
					$log->info("merged $merged and $file2");				
				}
				else {
					$log->info("rejecting $file2 from $merged");
				}		
			} else {
				$log->warn("profile alignment of $merged and $file2 failed")
			}
		}
		open my $outfh, '>>', $workdir . '/' . $outfile or die $!;
		print $outfh $merged, "\n";
		close $outfh;
	} @clusters;

	$log->info("DONE, results written to $outfile");
	return 1;
}


# Helper subroutines
sub _cluster {
	my ( $clusters, $hitset, $hitmap) = @_;
		
	# iterate over GIs in focal hit set
	for my $gi ( @{ $hitset } ) {
		
		# if focal GI has itself a set of hits, add those hits 
		# to current focal set and keep processing it
		if ( my $hits = delete $hitmap->{$gi} ) {
			_merge( $hitset, $hits );
			_cluster( $clusters, $hitset, $hitmap );
		}
	}		
	# when all GIs are processed, the grown hit set has become a cluster
	push @{ $clusters }, $hitset;	
	#$log->info("built single-linkage cluster ".scalar(@{ $clusters }));
}
	
sub _merge {
	my ( $set1, $set2 ) = @_;
	@$set1 = sort { $a <=> $b } keys %{ { map { $_ => 1 } ( @$set1, @$set2 ) } };
}

1;
