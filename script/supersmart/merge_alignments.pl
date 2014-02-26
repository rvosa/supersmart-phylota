#!/usr/bin/perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use List::Util 'reduce';
use Bio::SearchIO;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;

=head1 NAME

merge_alignments.pl - merges alignments by orthology

=head1 SYNOPSIS

 $ merge_alignments.pl --list=<list> --workdir=<dir> [--verbose] > <merged list>

=head1 DESCRIPTION

Given an input C<list> of alignment locations, assess their orthology (presently, on the
basis of all-vs-all BLAST searches) and profile align the orthologous clusters. The list
of so produced merged alignments is printed to STDOUT and the alignments themselves are
written to C<workdir>.

=cut

# process command line arguments
my $verbosity = WARN;
my ( $list, $workdir );
GetOptions(
	'verbose+'  => \$verbosity,
	'list=s'    => \$list,
	'workdir=s' => \$workdir,
);
my $dbname = File::Spec->catfile( $workdir, 'seeds.fa' );

# instantiate helper objects
my $service = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
my $mts     = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
my $config  = Bio::Phylo::PhyLoTA::Config->new;
my $log     = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# read list of alignments
my @alignments;
{
	$log->info("going to read file list $list");
	my %seen; # for some reason there might be duplicates
	open my $fh, '<', $list or die $!;
	while(<$fh>) {
		chomp;
		$seen{$_}++;
	}
	@alignments = keys %seen;
}

# write seed sequences
$log->info("going to write sequences to $dbname");
open my $fh, '>', $dbname or die $!;
for my $aln ( @alignments ) {
	if ( $aln =~ /(\d+)\.fa/ ) {
		my $gi = $1;
		my $seq = $service->find_seq($gi);
		print $fh '>', $gi, "\n", $seq->seq, "\n";
	}	
}

# make blast db
$log->info("going to make BLAST db for $dbname");
my @cmd = (
	$config->MAKEBLASTDB_BIN,
	'-in'     => $dbname,
	'-dbtype' => 'nucl',
	'2>'      => '/dev/null',
	'>'       => '/dev/null',		
); 
system("@cmd") and die "Problem making db: $?";

# run blast search
$log->info("going to run all vs all BLAST search on $dbname");
@cmd = (
	$config->BLASTN_BIN,
	'-query'  => $dbname,
	'-db'     => $dbname,
	'2>'      => '/dev/null',	
);
my $result = `@cmd`;

# process results, this is all-vs-all so many results with many hits. we will discard
# all hits where the overlapping region is smaller than 0.51 of the length of both
# query and hit
my $overlap = $config->MERGE_OVERLAP;
$log->info("going to process BLAST results");
my %hits;
my $q = 1;
open my $out, '<', \$result;
my $report = Bio::SearchIO->new( '-format' => 'blast', '-fh' => $out );
while ( my $result = $report->next_result ) {
	my $query = $result->query_name;
	my $q_l = length($service->find_seq($query)->seq);
	$hits{$query} = [];
	
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
				push @{ $hits{$query} }, $name;
				$log->debug("\thit: $name");
			}
			else {
				$log->debug("discarding hit $name for query $query");
			}
		}
	}
	
	# report progress
	my $template = 'found %i hits for query %s (%i nt) - processed %i / %i';
	my @args = ( scalar(@{$hits{$query}}), $query, $q_l, $q++, scalar(@alignments) );
	my $message  = sprintf $template, @args;
	$log->info($message);	
}

# make single linkage clusters
my $sets = [];
for my $gi ( keys %hits ) {
	if ( $hits{$gi} ) {
		$log->info("going to cluster from seed $gi");
		cluster( $sets, delete $hits{$gi}, \%hits );
	}
}

sub cluster {
	my ( $clusters, $hitset, $hitmap ) = @_;
	
	# iterate over GIs in focal hit set
	for my $gi ( @{ $hitset } ) {
	
		# if focal GI has itself a set of hits, add those hits 
		# to current focal set and keep processing it
		if ( my $hits = delete $hitmap->{$gi} ) {
			merge( $hitset, $hits );
			cluster( $clusters, $hitset, $hitmap );
		}
	}
	
	# when all GIs are processed, the grown hit set has become a cluster
	push @{ $clusters }, $hitset;	
	$log->info("built single-linkage cluster ".scalar(@{ $clusters }));
}

sub merge {
	my ( $set1, $set2 ) = @_;
	@$set1 = sort { $a <=> $b } keys %{ { map { $_ => 1 } ( @$set1, @$set2 ) } };
}

# now remove duplicates
@$sets = values %{ { map { join('|', sort { $a <=> $b } @{$_}) => $_ } @$sets } };

# align
my $i = 1;
my $total = scalar @{ $sets };
CLUSTER: for my $cluster ( @{ $sets } ) {
	$log->info("merging alignments in cluster $i / $total");
	
	# turn GIs into file names, check for singletons
	my @files = map { File::Spec->catfile( $workdir, $_ . '.fa' ) } @{ $cluster };
	if ( scalar(@files) == 1 ) {
		$log->info("singleton cluster $i: @files");
		print @files, "\n";
		$i++;
		next CLUSTER;
	}
	
	# the name of the file that contains the merger of this cluster
	my $merged = File::Spec->catfile( $workdir, "cluster${i}.fa" );
	
	# (the effectiveness of this procedure may depend on the input order). I was trying
	# to do this using List::Util::reduce, but it segfaults.
	$log->info("going to reduce @files");		
	my $file1;
	for my $i ( 0 .. ( $#files - 1 ) ) {
	
		# as long as merges are unsuccessful, we will continue to attempt
		# profile align the subsequence files against the first input file. TODO:
		# sort the input files by size so that if all mergers fail, we end up with
		# the largest input file
		if ( not $file1 ) {
			$file1 = $files[0];
		}
		my $file2 = $files[$i+1];
		
		# do the profile alignment
		$log->debug("attempting to merge $file1 and $file2");
		my $result = $service->profile_align_files($file1,$file2);
		
		# evaluate how this went
		my %fasta  = $mts->parse_fasta_string($result);
		if ( $mts->calc_mean_distance(%fasta) < $config->BACKBONE_MAX_DISTANCE ) {
			%fasta = $mts->dedup(%fasta);
			open my $fh, '>', $merged or die $!;
			for my $defline ( keys %fasta ) {
				print $fh '>', $defline, "\n", $fasta{$defline}, "\n";
			}
			$log->info("merged $file1 and $file2");
			
			# from now on we will continue to append to the merged file
			$file1 = $merged;
		}
		else {
			$log->info("rejecting $file2 from $file1");
		}		
	}
	print $merged || $file1, "\n";
	$i++;
}