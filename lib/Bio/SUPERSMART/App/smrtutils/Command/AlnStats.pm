package Bio::SUPERSMART::App::smrtutils::Command::AlnStats;

use strict;
use warnings;

use File::Spec;
use Data::Dumper;
use List::MoreUtils qw(uniq);

use Bio::Phylo::IO qw(parse);
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::PhyLoTA::Service::ParallelService;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: extract information about a set of alignments


=head1 NAME

AlnStats - Extract useful information given a list of alignment locations

=head1 SYNOPSYS

smrt-utils alnstats -a <file> [-o <file>] [-h] [-v] \ [-w <dir>] [-l <file>] [-y]

=head1 DESCRIPTION

The script produces two tables (tsv format). The first one,
giving sequence stats lists following alignment properties:

 - length of alignment
 - number of sequences in alignment
 - per cent gaps per sequence
 - number of insertions per alignment
 - number of deletions per alignment
 - average size of indels per alignment
 - number of invariant columns in alignment

The second one lists the alignments as columns and species as row. Species
present in alignment is indicated with '1' in the corresponding cell 

=cut

sub options {
	my ($self, $opt, $args) = @_;
	my $statstable_default = 'alignment-stats.tsv';
	my $markertable_default = 'marker-presence.tsv';
	return (
		['alignments|a=s', "list of alignment file locations", { arg => 'file', mandatory => 1 } ],
		["statstable|s=s", "name of output file with alignment stats, defaults to '$statstable_default'", {default => $statstable_default, arg => "file"}],
		["markertable|m=s", "name of output file with marker presences, defaults to '$markertable_default'", {default => $markertable_default, arg => "file"}],
	    );
}

sub validate {
	my ($self, $opt, $args) = @_;
	my $file = $opt->alignments;
	$self->usage_error("file $file does not exist") unless (-e $file);
	$self->usage_error("file $file is empty") unless (-s $file);
}

sub run {
	my ($self, $opt, $args) = @_;

	# get all alignment files
	open my $fh, '<', $opt->alignments or die $!;
	my @alnfiles = <$fh>;
	chomp @alnfiles;
	close $fh;

	# collect statistics for each matrix
	my @aln_stats;# = pmap {
#		$self->_get_aln_stats(@_);
#	} @alnfiles;

	my %marker_stats = $self->_get_marker_presence(@alnfiles);
	my @all_species = uniq map {keys(%{$_})} values( %marker_stats );
	
	
	open my $mfh, '>', $opt->markertable or die $!;
	print $mfh "species\t";
	print $mfh join("\t", @alnfiles);
	print $mfh "\n";
	for my $sp( @all_species ) {
		print $mfh $sp;
		for my $aln( @alnfiles ) {
			print $mfh "\t";
			my %h = %{$marker_stats{$aln}};
			print $mfh $h{$sp} ? 1: 0;

		}
		print $mfh "\n";
	}
   
	# write aln stats to table
	my @header = sort { lc ($a) cmp lc ($b) } keys %{$aln_stats[0]};
	open my $outfh, '>', $opt->statstable or die $!;
	print $outfh join("\t", @header) . "\n";
	for my $s ( @aln_stats ) {
		my %h = %{$s};
		print $outfh join("\t", @h{@header}) . "\n";		
	}


	$self->logger->info("DONE. Stats written to  " . $opt->statstable .  " marker presences written to " . $opt->markertable);
}

sub _get_marker_presence {
	my ($self, @alns) = @_;

	my %stats;
	for my $aln ( @alns ) {
		
		$stats{$aln} = {};
		
		# parse matrix
		my $project = parse(
			'-format'     => 'fasta',
			'-type'       => 'dna',
			'-file'     => $aln,
			'-as_project' => 1,
			);
		my ($matrix) = @{ $project->get_items(_MATRIX_) };
		
		for my $m( @{$matrix->get_entities} ) {
			my $species;
			if ( $m->get_name=~m/taxon\|(\d+)\//) {
				$species = $1;
			}
			else {
				$species = $m->get_name;
			}
			$stats{$aln}->{$species} = 1;
		}
	}
	return %stats;
}

# given a matrix, get properties of the alignment (gap count, indel count and sizes, etc)
sub _get_aln_stats  {
	my ($self, $aln) = @_;

	my %stats;
	$stats{'file'} = $aln;

	# parse matrix
	my $project = parse(
			'-format'     => 'fasta',
			'-type'       => 'dna',
			'-file'     => $aln,
			'-as_project' => 1,
			);
	my ($matrix) = @{ $project->get_items(_MATRIX_) };

	# number of seqs (taxa) and characters
	$stats{'nchar'} = $matrix->get_nchar;
	$stats{'ntax'} = $matrix->get_ntax;

	# get insertions and deletions
	my %deletions  = %{ $matrix->calc_indel_sizes( '-trim' => 1 ) };
	my %insertions = %{ $matrix->calc_indel_sizes( '-trim' => 1, '-insertions' => 1 ) };

	my @ins_sizes = keys %insertions;

	# count total number of insertions and deletions in alignment
	my $ins_count = 0;
	$ins_count += $_ for values(%insertions);
	my $del_count = 0;
	$del_count += $_ for values(%deletions);
	$stats{'ins_count'} = $ins_count;
	$stats{'del_count'} = $del_count;

	# average size of indels in alignment
	my $ins_avg_size = 0;
	$ins_avg_size += $_ for keys(%insertions);
	$ins_avg_size /= scalar(keys(%insertions)) if scalar(keys(%insertions));
	my $del_avg_size = 0;
	$del_avg_size += $_ for keys(%deletions);
	$del_avg_size /= scalar(keys(%deletions)) if scalar(keys(%deletions));
	$stats{'ins_avg_size'} = $ins_avg_size;
	$stats{'del_avg_size'} = $del_avg_size;

	# count gaps in all sequences
	my $gapcount = 0;
	my $avg_gapfreqs = 0;
	for my $dat (@{$matrix->get_entities}){
		my %freqs = %{ $matrix->calc_state_frequencies('-gap'=>1) };
		$avg_gapfreqs += $freqs{'-'} if $freqs{'-'};
		my $ch = $dat->get_char;
		my $cnt = $ch =~ tr/\-//;
		$gapcount += $cnt;
	}

	# average number of gaps per sequence
	$gapcount /= scalar(@{$matrix->get_entities});
	$stats{'gaps_per_seq'} = $gapcount;

	# average gap frequencies
	$avg_gapfreqs /= scalar(@{$matrix->get_entities});
	$stats{'gap_freq'} = $avg_gapfreqs;;

	# proportion of invariant sites
	my $prop_invar = $matrix->calc_prop_invar;
	$stats{'prop_invar'} = $prop_invar;

	return \%stats;
}

1;
