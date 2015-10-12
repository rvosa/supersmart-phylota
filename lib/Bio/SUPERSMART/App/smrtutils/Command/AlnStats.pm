package Bio::SUPERSMART::App::smrtutils::Command::AlnStats;

use strict;
use warnings;

use File::Spec;
use Data::Dumper;

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

The script produces a table (tsv format) listing following alignment properties:

 - length of alignment
 - number of sequences in alignment
 - per cent gaps per sequence
 - number of insertions per alignment
 - number of deletions per alignment
 - average size of indels per alignment
 - number of invariant columns in alignment

=cut

sub options {
	my ($self, $opt, $args) = @_;
	my $outfile_default = 'alignment-stats.tsv';
	return (
		['alignments|a=s', "list of alignment file locations", { arg => 'file', mandatory => 1 } ],
		["outfile|o=s", "name of output file with alignment stats, defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],
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
	my @stats = pmap {
		$self->_get_aln_stats(@_);
	} @alnfiles;

	# write stats to table
	my @header = sort { lc ($a) cmp lc ($b) } keys %{$stats[0]};
	open my $outfh, '>', $self->outfile or die $!;
	print $outfh join("\t", @header) . "\n";
	for my $s ( @stats ) {
		my %h = %{$s};
		print $outfh join("\t", @h{@header}) . "\n";		
	}


	$self->logger->info("DONE. Alignment stats printed to " . $self->outfile);
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
