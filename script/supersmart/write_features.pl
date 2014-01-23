#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;

=begin comment

We need to store the following:

| primary_tag  | varchar(12)          | YES  |     | NULL    |                |
| gi           | bigint(20) unsigned  | YES  | MUL | NULL    |                |
| gi_feat      | bigint(20) unsigned  | YES  | MUL | NULL    |                |
| ti           | bigint(20) unsigned  | YES  | MUL | NULL    |                |
| acc          | varchar(12)          | YES  |     | NULL    |                |
| acc_vers     | smallint(5) unsigned | YES  |     | NULL    |                |
| length       | bigint(20) unsigned  | YES  |     | NULL    |                |
| codon_start  | tinyint(3) unsigned  | YES  |     | NULL    |                |
| transl_table | tinyint(3) unsigned  | YES  |     | NULL    |                |
| gene         | text                 | YES  |     | NULL    |                |
| product      | text                 | YES  |     | NULL    |                |
| seq          | mediumtext           | YES  |     | NULL    |                |

In addition to that there is an auto-incrementing primary key. If we just create
a big text file (tab-delimited, with these headers in this order) we can then
import that in MySQL.

=cut comment

# process command line arguments
my ( $verbosity, $expand, $infile ) = ( WARN, 'zcat' );
GetOptions(
	'verbose+'  => \$verbosity,
	'infile=s'  => \$infile,
	'expand=s'  => \$expand,
);

# instantiate helper objects
my $sg  = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# instantiate file reader
$log->info("going to read genbank records from '$expand $infile'");
open my $fh, '-|', "$expand $infile" or die $!;
my $reader = Bio::SeqIO->new(
	'-format' => 'genbank',
	'-fh'     => $fh,
);

# print header
my @columns = (
	'primary_tag',  # bioperl provides this
	'gi',           # this is the foreign key to the seqs table
	'gi_feat',      # the gi of the feature itself
	'ti',           # the taxon of the feature
	'acc',          # accession number
	'acc_vers',     # accession version number
	'length',       # length of the subsequence, if applicable
	'codon_start',  # reading frame start, if applicable
	'transl_table', # which translation table applies, if any
	'gene',         # the gene name, i.e. the important bit for us, e.g. rbcL
	'product',      # textual description of the gene
	'seq',          # the subsequence
);
print join("\t", @columns), "\n";

# start iterating over the records
while ( my $seq = $reader->next_seq ) {
	
	# first check to see if the sequence is in the seqs table
	my $gi = $seq->primary_id;
	if ( $sg->find_seq($gi) ) {
		$log->info("GI:$gi is in phylota");
				
		# iterate over all
		my $ti;
		FEATURE: for my $feat ( $seq->get_SeqFeatures ) {
			my $primary_tag = $feat->primary_tag;
			
			# some types of features we should skip as they're useless
			next FEATURE if $primary_tag eq 'misc_feature';
			next FEATURE if $primary_tag eq 'exon';
			next FEATURE if $primary_tag eq 'intron';
			next FEATURE if $primary_tag eq 'rRNA';
			next FEATURE if $primary_tag eq 'sig_peptide';
			next FEATURE if $primary_tag eq 'mat_peptide';
			next FEATURE if $primary_tag eq 'polyA_signal';
			next FEATURE if $primary_tag eq '5\'UTR';
			next FEATURE if $primary_tag eq 'gene';
			
			# this is going to be a pseudo object to write to file
			my %feature = ( primary_tag => $primary_tag, gi => $gi );
			$log->info("Primary tag: $feature{primary_tag}");			
			
			# iterate over tags			
			for my $tag ( $feat->get_all_tags ) {
				for my $value ( $feat->get_tag_values($tag) ) {
					
					# tag is a reference to another database record, i.e. the
					# GI of a subsequence, or a taxon
					if ( $tag eq 'db_xref' ) {
						if ( $value =~ /^taxon:(\d+)$/ ) {
							$feature{ti} = $1;
							$ti = $feature{ti};
						}
						elsif ( $value =~ /^GI:(\d+)$/ ) {
							$feature{gi_feat} = $1;
						}
					}
					
					# value is the protein translation
					elsif ( $tag eq 'translation' ) {
						$feature{seq} = $value;
						$feature{length} = length $value;
					}
					
					# value is the gene name
					elsif ( $tag eq 'gene' ) {
						$feature{gene} = $value;
					}
					
					# value is the textual description of the locus
					elsif ( $tag eq 'product' ) {
						$feature{product} = $value;
					}
					
					# value is the open reading frame of the gene
					elsif ( $tag eq 'codon_start' ) {
						$feature{codon_start} = $value;
					}
					
					# value is the accession number
					elsif ( $tag eq 'protein_id' && $value =~ /^(.+?)\.(\d+)$/ ) {
						my ( $acc, $ver ) = ( $1, $2 );
						$feature{acc} = $acc;
						$feature{acc_vers} = $ver;
					}
					
					$log->info("\t$tag => $value");
				}
			}
			
			# here we do the print that is re-directed to the output file for
			# import in mysql. The schema seems to want a record for each
			# primary tag, which means the records will be a bit sparse.
			# We therefore want to silence the 'unitialized' warnings.
			if ( $primary_tag eq 'CDS' ) {
				$feature{ti} = $ti;
				{
					no warnings 'uninitialized';
					print join("\t", @feature{@columns}), "\n";
					use warnings;
				}
			}
		}
	}
	
	# the sequence is not in the table, skip
	else {
		$log->warn("GI:$gi is not in phylota, skipping...");
	}
}
