#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
#use Parallel::MPI::Simple;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::DB::GenBank;
use Bio::Phylo::Util::Logger ':levels';


my ( $verbosity ) = ( WARN );
GetOptions( 'verbose+'  => \$verbosity);

# instantiate helper objects
my $sg  = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

my $db_obj = Bio::DB::GenBank->new;

# get all gi numbers from sequence table
my $schema = 'Bio::Phylo::PhyLoTA::DAO'->new;
my $result = $schema->resultset('Seq')->search_rs({}, {columns=>'gi'});

my $count = $result->count;
#print "Count = $count \n";

my @slices;

my @sl = $result->slice(0, 1000); #caution magic number
my @gis = map($_->gi, @sl);


# make equally sized chunks of lists of GI codes
#my $counter = 0;
#while($counter < $count){
#    push @slices, 
#}

# get sequence objects from NCBI server
my $stream = $db_obj->get_Stream_by_id(\@gis);

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

my $cnt = 0;
while(my $seq = $stream->next_seq){
    my $gi = $gis[$cnt];
    $cnt++;
    #print "GI: $gi \n";
    FEATURE: for my $feat ( $seq->get_SeqFeatures ) {		
	my $primary_tag = $feat->primary_tag;
	#print "Primary tag : $primary_tag \n";       

	my $ti;
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







#my $fac = Bio::DB::SoapEUtilities->new();

#my $res = $fac->esearch(
#    -db => 'nucleotide',
#    -ids=>\@gis,
#    -usehistory => 1 )->run( -no_parse => 1 );

##my $query = Bio::DB::Query::GenBank->new(-db => 'nucleotide',
##					 -ids => \@gis,
##					 -WebEnv => $res->webenv,
##					 -QueryKey => $res->query_key);

#use Bio::DB::SoapEUtilities;

##print "Getting seqs\n";
##my $stream = $db_obj->get_Stream_by_query($query);
##print "Seqs retrieved\n";




