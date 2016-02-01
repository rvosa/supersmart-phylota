package Bio::SUPERSMART::App::smrtutils::Command::AlnStats;

use strict;
use warnings;

use File::Spec;
use Data::Dumper;
use List::MoreUtils qw(uniq);

use Bio::Phylo::IO qw(parse);
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::SUPERSMART::Service::ParallelService;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: extract information about a set of alignments

=head1 NAME

AlnStats - Extract useful information given a list of alignment locations

=head1 SYNOPSYS

smrt-utils alnstats -a <file> [-o <file>] [-h] [-v] \ [-w <dir>] [-l <file>] [-y]

=head1 DESCRIPTION

The script produces a tables (tsv format)
giving sequence stats lists following alignment properties:

 - length of alignment
 - number of sequences in alignment
 - per cent gaps per sequence
 - number of insertions per alignment
 - number of deletions per alignment
 - average size of indels per alignment
 - number of invariant columns in alignment
 - NCBI identifiers for species in the table

=cut

sub options {
	my ($self, $opt, $args) = @_;
	my $statstable_default = 'alignment-stats.tsv';
	return (
		['alignments|a=s', "list of alignment file locations", { arg => 'file', mandatory => 1 } ],
		["statstable|s=s", "name of output file with alignment stats, defaults to '$statstable_default'", {default => $statstable_default, arg => "file"}],
	    );
}

sub validate {
	my ($self, $opt, $args) = @_;
	my $file = $opt->alignments;
	$self->usage_error("input file not given") unless ( $file);
	$self->usage_error("input file does not exist") unless (-e $file);
	$self->usage_error("file $file is empty") unless (-s $file);
}

sub run {
	my ($self, $opt, $args) = @_;

	# get alignment files
	open my $fh, '<', $opt->alignments or die $!;
	my @alnfiles = <$fh>;
	chomp @alnfiles;
	close $fh;

	# collect statistics for each matrix
	my @aln_stats = pmap {
		$self->_get_aln_stats(@_);
	} @alnfiles;

	# write table with statistics on all alignments
	my %h = %{$aln_stats[0]};
	my @header = sort { lc ($a) cmp lc ($b) } keys (%h);
	open my $outfh, '>', $opt->statstable or die $!;
	print $outfh join("\t", @header) . "\n";
	for my $s ( @aln_stats ) {
		my %hash = %{$s};
		print $outfh join("\t", @hash{@header}) . "\n";
	}

	$self->logger->info("DONE. Stats written to  " . $opt->statstable);
}

# given a matrix, get properties of the alignment (gap count, indel count and sizes, etc)
sub _get_aln_stats  {
	my ($self, $file) = @_;

	# parse matrix
	my $project = parse(
		'-format'     => 'fasta',
		'-type'       => 'dna',
		'-file'     => $file,
		'-as_project' => 1,
		);

	my ($matrix) = @{ $project->get_items(_MATRIX_) };

	my %stats;
	$stats{'file'} = $file;

	# number of seqs (taxa) and characters
	$stats{'nchar'} = $matrix->get_nchar;
	$stats{'ntax'} = $matrix->get_ntax;

	# get insertions and deletions
	my %deletions  = %{ $matrix->calc_indel_sizes( '-trim' => 0 ) };
	my %insertions = %{ $matrix->calc_indel_sizes( '-trim' => 0, '-insertions' => 1 ) };

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

	my @species;
	# count gaps in all sequences
	my $gapcount = 0;
	my $avg_gapfreqs = 0;
	for my $dat (@{$matrix->get_entities}){
		my %freqs = %{ $matrix->calc_state_frequencies('-gap'=>1) };
		$avg_gapfreqs += $freqs{'-'} if $freqs{'-'};
		my $ch = $dat->get_char;
		my $cnt = $ch =~ tr/\-//;
		$gapcount += $cnt;

		# get species id
		my $spec;
		if ( $dat->get_name=~m/taxon\|(\d+)/) {
			$spec = $1;
		}
		else {
			$spec = $dat->get_name;
		}
		push @species, $spec;
	}
	$stats{'species'} = join(',', @species);

	# average number of gaps per sequence
	$gapcount /= scalar(@{$matrix->get_entities});
	$stats{'gaps_per_seq'} = $gapcount;

	# average gap frequencies
	$avg_gapfreqs /= scalar(@{$matrix->get_entities});
	$stats{'gap_freq'} = $avg_gapfreqs;;

	# proportion of invariant sites
	my $prop_invar = $matrix->calc_prop_invar;
	$stats{'prop_invar'} = $prop_invar;

	# get number of identical sequence pairs
	my $ident_pairs = $self->_num_identical_seqs($matrix);
	$stats{'ident_pairs'} = $ident_pairs;

	# get average distance of seqs in alignment
	my $avg_dist = $self->_avg_dist($file);
	$stats{'avg_dist'} = $avg_dist;

	return \%stats;
}

sub _avg_dist {
	my ($self, $filename) = @_;
	my $mt    = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
	open my $fh, '<', $filename;
	read $fh, my $string, -s $fh;
	close $fh;
	my $dist = $mt->calc_mean_distance($string);
	return $dist;
}

sub _num_identical_seqs {
	my ($self, $matrix) = @_;

	my $mts    = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;

	my @rows = @{$matrix->get_entities};

	my $identical = 0;
	for (my $i=0; $i<=$#rows; $i++) {
		my $row_i = $rows[$i];
		for (my $j=$i+1; $j<=$#rows; $j++) {
			my $row_j = $rows[$j];

			if ($row_i->get_char eq $row_j->get_char) {

				my $name_i = $row_i->get_name;
				my $name_j = $row_j->get_name;
				if ( $name_i=~/.*taxon\|([0-9]+)/) {
					$name_i =$1;
				}
				if ( $name_j=~/.*taxon\|([0-9]+)/) {
					$name_j =$1;
				}

				$name_i = $mts->find_node($name_i)->taxon_name;
				$name_j = $mts->find_node($name_j)->taxon_name;

				$identical++;
				$self->logger->info("seqs " . $name_i . " and " . $name_j  . " are identical !");
			}
		}
	}
	return $identical;
}



1;
