# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Service::SequenceGetter; # maybe this should be SequenceService
use strict;
use warnings;
use Moose;
use Data::Dumper;
use File::Temp 'tempfile';
use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::Phylo::Factory;
use Bio::Phylo::Util::Exceptions 'throw';
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Phylo::IO 'parse_matrix';

extends 'Bio::Phylo::PhyLoTA::Service';

my $fac = Bio::Phylo::Factory->new;

=head1 NAME

Bio::Phylo::PhyLoTA::Service::SequenceGetter - Sequence Getter

=head1 DESCRIPTION

This package interacts with sequence data and clusters. XXX: this should probably
become something like SequenceService, or AlignmentService.

=head1 METHODS

=over

=item store_genbank_sequences

This stores genbank sequences from a file.

=cut

sub store_genbank_sequences {
    my ( $self, $file ) = @_;
    
    # also allow gzipped files, pipe from gunzip -c
    my $command = ( $file =~ /\.gz$/ ) ? "gunzip -c $file |" : "< $file";
    
    # instantiate bioperl sequence reader
    my $reader  = Bio::SeqIO->new(
        '-format' => 'genbank',
        '-file'   => $command,
    );
    
    # iterate over sequences in file
    while( my $seq = $reader->next_seq ) {
        $self->store_sequence($seq);
    }
}

=item store_sequence

This stores a bioperl sequence object, if optional second argument is true we
won't store the raw sequence string (because it would be too long). Returns
the newly created dao Seq object.

=cut

sub store_sequence {
    my ( $self, $bioperl_seq, $no_raw ) = @_;
    my @dates = $bioperl_seq->get_dates();
    return $self->schema->resultset('Seq')->create( {
        'acc_date' => _formatDate( $dates[-1] ),
        'gbrel'    => $self->config->currentGBRelease,
        'seq'      => $no_raw ? undef : $bioperl_seq->seq,        
        'acc'      => $bioperl_seq->accession_number,
        'acc_vers' => $bioperl_seq->seq_version,
        'def'      => $bioperl_seq->desc,
        'division' => $bioperl_seq->division,
        'gi'       => $bioperl_seq->primary_id,
        'length'   => $bioperl_seq->length,
        'ti'       => $bioperl_seq->species->ncbi_taxid
    });
}

=item store_feature

This populate the CDS and RNA features, taking care with remotely accessioned features.
NB! Bioperl feature->spliced_seq will just return a guess at the length of the sequence, padded with 'N's
when the acc number is remote. This is often a bad guess because it is based on the presumption that the
ENTIRE feature is remote, when often just a piece of the feature is remote. Go ahead, look at the code...
if( !defined $called_seq ) {
	$seqstr .= 'N' x $self->length;  ...here the length is for the feature's location, not the features sublocation
	next; ...so for something like join(BC123.1:1-100, 12-200,10000-10100) it might be 10100 minus 12.
DO NOT USE feat->length for split sequences at all! unless you want the length of the whole region from min to max.

=cut

sub store_feature {
    my ( $self, $daoseq, $feat, $no_raw ) = @_;

    # parameters for the Feature create method below
    my %params = (
        'primary_tag'   => $feat->primary_tag,
		'gi'            => $daoseq->gi,
		'ti'            => $daoseq->ti,
        'range'         => $feat->location->to_FTstring(),
        
        # these we will populate out of the tags we read below
        'gene'          => undef, # gene symbol
        'transl_table'  => undef, # translation table
        'codon_start'   => undef, # reading frame
        'product'       => undef, # gene product
        'gi_feat'       => undef, # db_xref        
        'acc'           => undef, # accession number, the part before the version
        'acc_vers'      => undef, # accession number version, i.e. a number after the dot
        'seq'           => undef, # the raw sequence, may remain undef if remote
        'length'        => undef, # length of the raw sequence
    );    
    
    # if we don't trap for remote sequences, bad things happen, see above.
    if ( ! grep { $_->is_remote } $feat->location->each_Location() && ! $no_raw ) {
        $params{'seq'} = $feat->spliced_seq()->seq();
        $params{'length'} = length $params{'seq'};
    }
    
    # every feature has a number of key value pairs associated with them
    for my $tag ( $feat->get_all_tags() ) {
        my $value = join( ' ', $feat->get_tag_values($tag) );
        $params{$tag} = $value if exists $params{$tag};
        
        # only keep the numerical part of db_xref
        if ( $tag eq 'db_xref' && $value =~ /(\d+)/ ) {
            $params{'gi_feat'} = $1;
        }
        
        # split protein_id as an accession number
        if ( $tag eq 'protein_id' && $value =~ /^([^.]+)\.(.*)$/ ) {
            ( $params{'acc'}, $params{'acc_vers'} ) = ( $1, $2 );
        }
    }
    
    # some CDS features are not taken seriously, not translated, etc., so skip
    if ( $feat->primary_tag =~ /CDS/ ) {        
        return if not defined $params{'acc'};
    }

    return $self->schema->resultset('Feature')->create( \%params );
}

=item _formatDate

This format dates...

=cut

sub _formatDate {
    my %monthH = (
        JAN => 1,
        FEB => 2,
        MAR => 3,
        APR => 4,
        MAY => 5,
        JUN => 6,
        JUL => 7,
        AUG => 8,
        SEP => 9,
        OCT => 10,
        NOV => 11,
        DEC => 12
    );
    my ( $day, $month, $year ) = split '\-', $_[0];
    return "$year-$monthH{$month}-$day";
}

=item get_aa_for_sequence

This fetches the amino acid sequence (if any) for the provided nucleotide GI.

=cut

sub get_aa_for_sequence {
	my ($self,$gi) = @_;
	my $gb = Bio::DB::GenBank->new;
	
	# see Bio::DB::WebDBSeqI, this needs to be set to string because otherwise
	# we fork a new process, and these forked children aren't being cleaned up so
	# they accumulate over time
	$gb->retrieval_type('io_string');
	
	if ( my $seq = $gb->get_Seq_by_gi($gi) ) {
		
		# iterate over sequence features
		for my $feat ( $seq->get_SeqFeatures ) {
			
			# check to see if it's a protein coding sequence with an AA translation
			if ( $feat->primary_tag eq 'CDS' and $feat->has_tag('translation') ) {
				
				# fetch and pring the translation
				my ($protseq) = $feat->get_tag_values('translation');
				return $protseq;
			}
		}		
	}
	return;
}

=item run_blast_search

This runs a local BLAST search. We initially assume we run the search against
the InParanoid data base, which is defined in the config object, and that
we run a protein blast, i.e. 'blastp'. The method needs a -seq => $seq argument,
which can either be a Seq object, or a raw string such as is produced by
get_aa_for_sequence. What we get back is the following: every BLAST hit will
have a sequence ID for which multiple InParanoid records exist, because that
sequence has been blasted against 99 other genomes, so the return value is
an array of result sets.

If we are comparing two root clusters, and for each of those we get the seed_gi,
we need to then get the amino acid translation for that, BLAST it against the
database, and then in the first result set, there has to be a combination of
results, i.e. InParanoid objects, for which $r1->is_orthologous($r2) is true.

=cut

sub run_blast_search {
	my ( $self, %args ) = @_;
	my $log = $self->logger;
	
	# at a minimum, we must have a -seq argument
	if ( not $args{'-seq'} ) {
		throw 'BadArgs' => "need -seq argument";
	}
	
	# need the seq as a raw string, not an object, so check to see if it's
	# a reference. If it is, call the ->seq method to get the raw data.
	# otherwise, just copy the value (i.e. the seq string).
	my $seq = ref $args{'-seq'} ? $args{'-seq'}->seq : $args{'-seq'};
	
	# create a temporary file with the sequence
	my ( $fh, $filename ) = tempfile();
	print $fh ">DUMMY\n$seq\n";
	close $fh;
	
	# run blast+
	my $report = Bio::Tools::Run::StandAloneBlastPlus->new(
		'-db_name' => $args{'-database'} || $self->config->INPARANOID_SEQ_FILE,
 	)->run(
 		'-method' => args{'-program'}  || 'blastp',
 		'-query'  => $filename,
	);
	
	# these will be InParanoid database objects	
	my @hits;	
	while ( my $result = $report->next_result() ) {
		$log->debug("iterating over result $result");
		
		# there should be just one result with several hits
		while ( my $hit = $result->next_hit() ) {
			if ( my $name = $hit->name ) {
				
				# there are odd suffixes in the file we downloaded from
				# inparanoid
				$name =~ s/_+spec_id_\d+$//;
				$log->debug("iterating over hit $hit");
				$log->debug("hit name is: $name");
				push @hits, $self->search_inparanoid({ 'protid' => $name });
			}
		}
	}
	$log->debug("found ".scalar(@hits)." hits");
	return @hits;	
}

=item get_orthologs_for_protein_id

Fetches all orthologous InParanoid protein IDs for a given ID

=cut

sub get_orthologs_for_protein_id {
    my ($self,$protid) = @_;
    my $log = $self->logger;
    my @orthologs;
    
    # this fetches all instances where $protid is a member in a pairwise comparison
    my @records = $self->search_inparanoid({ 'protid' => $protid })->all;
    $log->info("found ".scalar(@records)." records for $protid");
    my $i = 1;
    for my $record ( @records ) {
	$log->info("processing hit $i for protid $protid");
	
	# together, these two uniquely identify the pairwise comparison
	my $id = $record->id;
	my $guid = $record->guid;
	$log->info("related records will have ID $id and GUID $guid");
	
	# this gets all participants in the pairwise comparison (i.e. also our query $protid)
	my @results = $self->search_inparanoid({ 'id' => $id, 'guid' => $guid, 'bootstrap' => '100%' })->all;
	$log->info("found ".scalar(@results)." related records with ID $id and GUID $guid");
	
	# here we filter out the results where the member *is* $protid, keeping the other member
	push @orthologs, grep { $_ ne $protid } map { $_->protid } @results;
	$i++;
    }
    return @orthologs;
}

=item get_protid_for_seed_gi

Takes a seed GI, attempts to locate a protein translation for it,
BLASTs the protein against InParanoid and returns the protein ID
of the best hit

=cut

sub get_protid_for_seed_gi {
    my ( $self, $seed_gi ) = @_;
    my $log = $self->logger;
    $log->debug("going to analyse seed GI $seed_gi");
    my $protid;
	
    # do protein translation
    eval {
        if ( my $aa = $self->get_aa_for_sequence($seed_gi) ) {
            $log->debug("seed GI $seed_gi has protein translation: $aa");
            
            # run blast
            if ( my @hits = $self->run_blast_search( '-seq' => $aa ) ) {
                    
                # get best inparanoid hit
                if ( $hits[0] and $hits[0]->count > 0 ) {
                    $log->debug("seed GI $seed_gi has BLAST hits");
                    
                    # fetch and store protein id
                    $protid = $hits[0]->next->protid;                    
                    $log->debug("seed GI $seed_gi has best hit $protid");

                }
            }
        }
    };
    if ( $@ ) {
        $log->warn("couldn't fetch AA for sequence $seed_gi: $@");
    }
    return $protid;
}

=item get_largest_cluster_for_sequence

This fetches the largest containing cluster for the seed sequence. This is a higher
taxon, such as an order, for example.

=cut

sub get_largest_cluster_for_sequence{
    my($self,$gi)=@_;
    
    # get the logger from Service.pm, the parent class
    my $log=$self->logger;
    
    # search CiGi table for all subtree clusters that include provided gi
    my $cigis = $self->schema->resultset('CiGi')->search({ gi => $gi, cl_type => "subtree" });
    
    # clusterid for most inclusive cluster
    my $biggestcluster;
    
    # size of the most inclusive cluster
    my $clustersize=0;
    
    # root_taxon of the most inclusive cluster
    my $taxonid;
    
    # iterate over search results
    while(my $c=$cigis->next){
	$log->info($c->clustid," ",$c->ti,"\n");
	
	# search cluster table for all clusters with clustid and ti from cigi result
	my $clusters=$self->schema->resultset('Cluster')->search({ ci => $c->clustid, ti_root => $c->ti});
	
	# iterate over search results
	while (my $cluster=$clusters->next){
	    
	    # looking for most inclusive cluster with largest n_ti
	    if ($cluster->n_ti > $clustersize ){
		$clustersize=$cluster->n_ti;
		$biggestcluster=$cluster->ci;
		$taxonid=$cluster->ti_root;
	    }
	    $log->info("CLUSTER: ",$cluster->pi," ",$cluster->n_ti,"\n");
	}
    }
    $log->info("biggestcluster: ",$biggestcluster," ",$taxonid->ti,"\n");
    
    # search CiGi table for sequences with most inclusive cluster
    return $self->get_sequences_for_cluster_object({ cl_type => 'subtree', ti => $taxonid->ti, clustid => $biggestcluster});
}

=item get_smallest_cluster_for_sequence

This method fetches the smallest containing cluster for the seed sequence, typically
these are the results of all-vs-all blasting within a fairly low level taxon such as
a family or a large genus (according to the NCBI taxonomy).

=cut

sub get_smallest_cluster_for_sequence{
    my($self,$gi)=@_;
    
    # get the logger from Service.pm, the parent class
    my $log=$self->logger;
    
    # search CiGi table for all subtree clusters that include provided gi
    my $cigis = $self->schema->resultset('CiGi')->search({ gi => $gi, cl_type => "subtree" });
    
    # clusterid for least inclusive cluster
    my $smallestcluster;
    
    # size of the least inclusive cluster
    my $clustersize=undef;
    
    # root_taxon of the least inclusive cluster
    my $taxonid;
    
    # iterate over search results
    while(my $c=$cigis->next){
	$log->info($c->clustid," ",$c->ti,"\n");
	
	# search cluster table for all clusters with clustid and ti from cigi result
	my $clusters=$self->schema->resultset('Cluster')->search({ ci => $c->clustid, ti_root => $c->ti});
	
	# iterate over search results
	CLUSTER: while (my $cluster=$clusters->next){
	    
	    if (not defined $clustersize){
		$clustersize=$cluster->n_ti;
		$smallestcluster=$cluster->ci;
		$taxonid=$cluster->ti_root;
		next CLUSTER;
	    }
	    
	    # looking for least inclusive cluster with largest n_ti
	    if ($cluster->n_ti < $clustersize ){
		$clustersize=$cluster->n_ti;
		$smallestcluster=$cluster->ci;
		$taxonid=$cluster->ti_root;
	    }
	    $log->info("CLUSTER: ",$cluster->pi," ",$cluster->n_ti,"\n");
	}
    }
    $log->info("smallestcluster: ",$smallestcluster," ",$taxonid->ti,"\n");
    
    # search CiGi table for sequences with most inclusive cluster
    return $self->get_sequences_for_cluster_object({ cl_type => 'subtree', ti => $taxonid->ti, clustid => $smallestcluster});
}

=item get_sequences_for_cluster_object

This fetches the sequences for the cluster object. The $cluster_object is a hash
reference with the following keys:
    # - cl_type => either node or subtree
    # - ti      => the taxon id for the root taxon
    # - clustid => the cluster id, which is NOT a primary key

=cut

sub get_sequences_for_cluster_object {
    # $cluster_object is a hash reference with the following keys:
    # - cl_type => either node or subtree
    # - ti      => the taxon id for the root taxon
    # - clustid => the cluster id, which is NOT a primary key
    my ($self,$cluster_object) = @_;
    
    # search CiGi table to fetch all GIs within this cluster
    my $gis = $self->search_ci_gi($cluster_object);

    # this will hold the resulting sequences
    my @sequences;
    
    # iterate over search results
    while(my $gi = $gis->next){
	
	# look up sequence by it's unique id
	my $seq = $self->find_seq($gi->gi);
	
	# add sequence to results
	push @sequences, $seq;
    }
    
    # return results
    return @sequences;    
}

=item get_smallest_cluster_object_for_sequence

This fetches the smallest, least inclusive cluster object for provided gi.

=cut

sub get_smallest_cluster_object_for_sequence{
    my($self,$gi)=@_;
    
    # get the logger from Service.pm, the parent class
    my $log=$self->logger;
    
    # search CiGi table for all subtree clusters that include provided gi
    my $cigis = $self->search_ci_gi({ gi => $gi, cl_type => "subtree" });
    
    # clusterid for least inclusive cluster
    my $smallestcluster;
    
    # size of the least inclusive cluster
    my $clustersize=undef;
    
    # root_taxon of the least inclusive cluster
    my $taxonid;
    
    # iterate over search results
    while(my $c=$cigis->next){
	$log->info($c->clustid," ",$c->ti,"\n");
	
	# search cluster table for all clusters with clustid and ti from cigi result
	my $clusters=$self->schema->resultset('Cluster')->search({ ci => $c->clustid, ti_root => $c->ti});
	
	# iterate over search results
	CLUSTER: while (my $cluster = $clusters->next){
	    
	    if (not defined $clustersize){
		$clustersize = $cluster->n_ti;
		$smallestcluster = $cluster->ci;
		$taxonid = $cluster->ti_root;
		next CLUSTER;
	    }
	    
	    # looking for least inclusive cluster with largest n_ti
	    if ($cluster->n_ti < $clustersize ){
		$clustersize = $cluster->n_ti;
		$smallestcluster = $cluster->ci;
		$taxonid = $cluster->ti_root;
	    }
	    $log->info("CLUSTER: ",$cluster->pi," ",$cluster->n_ti,"\n");
	}
    }
    $log->info("smallestcluster: ",$smallestcluster," ",$taxonid->ti,"\n");
    
    # input to search CiGi table for sequences within cluster
    return { cl_type => 'subtree', ti => $taxonid->ti, clustid => $smallestcluster };
}

=item get_parent_cluster_object

This fetches the parent clusters of provided clusterid.

=cut

sub get_parent_cluster_object {
    my ($self,$cluster) = @_;
    my $ci = $cluster->{clustid} || $cluster->{ci};
    my $ti = $cluster->{ti} || $cluster->{ti_root};
    my $node = $self->schema->resultset('Node')->find($ti);
    my $parent_id = $node->ti_anc;
    return { cl_type => 'subtree', ti => $parent_id, clustid => $ci };
}

=item get_child_cluster_objects

This fetches the child clusters for provided clusterid.

=cut

sub get_child_cluster_objects {
    my ($self,$cluster) = @_;
    my $ci = $cluster->{clustid} || $cluster->{ci};
    my $ti = $cluster->{ti} || $cluster->{ti_root};
    my $nodes = $self->schema->resultset('Node')->search({ ti_anc => $ti });
    my @result;
    while(my $node = $nodes->next) {
	push @result, { cl_type => 'subtree', ti => $node->ti, clustid => $ci };
    }
    return @result;
}    

=item compute_median_seq_length

This computes the most occuring sequence length, i.e. the median length, within the cluster.

=cut

sub compute_median_seq_length {
    my ($self,@seq) = @_;
    
    # count number of occurrences of each sequence length
    my %occurrences;
    for my $seq (@seq) {
	$occurrences{ $seq->length }++;
    }
    
    # similar to $clustersize (line 141) and $biggestcluster (line 138)
    my ($n_occurrences,$most_seen_length) = (0,0);
    
    # iterate over all lengths
    for my $length ( keys %occurrences ) {
	
	# similar to line 156
	if ( $occurrences{$length} > $n_occurrences ) {
	    $most_seen_length = $length;
	    $n_occurrences = $occurrences{$length};
	}
    }
    
    # e.g. for Cytochrome B this should return 1140
    return $most_seen_length;
}

=item get_aligned_locus_indices

Given a GI, a locus name and an aligned sequence string, returns an array of zero-based 
indices that specify columns that occur inside a named locus.

=cut

sub get_aligned_locus_indices {
	my ( $self, $gi, $locus, $seq ) = @_;
	my $feat   = $self->single_feature({ gene => $locus, gi => $gi });
	my $start  = $feat->codon_start;
	my $length = $feat->length;
	my $end    = $start + $length;
	my $nucs   = 0;
	my @seq    = split //, $seq;
        my @indices;
	for my $i ( 0 .. $#seq ) {
		$nucs++ if $seq[$i] ne '-';
		push @indices, $i if $nucs >= $start && $nucs <= $end;
	}
	return @indices;
}

=item filter_seq_set

This method filters out duplicate sequences per taxon, and prefers to keep a
sequence of the median length accross the whole set, otherwise keep the shortest
of the ones that are longer than that or the longest of the ones that are
shorter.

=cut

# TODO: make it so that all else being equal we prefer the sequence with the lowest number of NNNN

sub filter_seq_set {
    my ($self,@seq) = @_;
    
    # compute accross the whole set
    my $median_length = $self->compute_median_seq_length(@seq);
    
    # this is what we return, sequence objects
    my @filtered_seqs;
    
    # multidimensional hash, top level keys are taxon ids, second level keys are raw sequences, values are sequence
    # objects
    my %sequences_for_taxon;
    
    # iterate over all sequences
    for my $seq(@seq) {
	
		# this fetches the NCBI taxon id
		my $ti = $seq->ti;
		
		# create an anonymous hash first time we see this taxon id
		if(not exists $sequences_for_taxon{$ti}){
			$sequences_for_taxon{$ti}={};
		}
		
		# assign the sequence object to raw sequence as key, this filters out duplicate sequences
		$sequences_for_taxon{$ti}->{$seq->seq}=$seq;
    }
    
    # iterate over taxon ids
    for my $ti(keys %sequences_for_taxon){
	
		# contains raw sequence strings after filtering by length
		my $raw_seq;
		
		# contains raw sequence strings before filtering by length
		my @taxon_seqs = keys %{ $sequences_for_taxon{$ti} };
		
		# only keep sequences of median length
		my @median_length_seqs = grep { length($_) == $median_length } @taxon_seqs;
		
		# only keep sequences larger than median length
		my @longer_seqs = sort { length($a) <=> length($b) } grep { length($_) > $median_length } @taxon_seqs;
						
		# tests if there were sequences of median length
		if ( @median_length_seqs ) {
			$raw_seq = $median_length_seqs[0];
		}
		
		# tests if there were sequences larger than median length
		elsif ( @longer_seqs ) {
			$raw_seq = $longer_seqs[0];
		}
		
		# there were only shorter sequences
		else {
			my @shorter_seqs = sort { length($b) <=> length($a) } @taxon_seqs;
			$raw_seq = $shorter_seqs[0];
		}
		
		# map raw sequences back to sequence objects
		push @filtered_seqs,$sequences_for_taxon{$ti}->{$raw_seq};
    }
    
    # return results
    return @filtered_seqs;
}

=item align_sequences

This method align sequences using the bioperl wrapper for muscle.

=cut

sub align_sequences{
    my ($self,@seq)=@_;
    
    # here we convert sequence objects from the database, i.e. Bio::Phylo::PhyLoTA::DAO::Result::Seq
    # objects (which are not compatible with bioperl) into Bio::Phylo::Matrices::Datum object, which
    # ARE compatible, so that we can pass those into the muscle wrapper
    my @convertedseqs;
    for my $seq(@seq){
		my $converted=$fac->create_datum(
			'-type' => 'dna',
			'-name' => $seq->gi,
			'-char' => $seq->seq,
		);
		push @convertedseqs,$converted;
    }
    
    # instantiate using factory, to address milestone issue #4
	my $aligner = $self->_make_aligner;
    my $alignment = $aligner->align(\@convertedseqs);
    return $alignment;
}

sub _muscle {
	my $muscle = shift;
	my @in = @Bio::Tools::Run::Alignment::Muscle::MUSCLE_SWITCHES;
	my @out = grep { $_ ne 'profile' } @in;
	@Bio::Tools::Run::Alignment::Muscle::MUSCLE_SWITCHES = @out;
	$muscle->quiet(1);
}

sub _quiet { shift->quiet(1) }

sub _verbose { shift->verbose(-1) }

sub _make_aligner {
	my $self = shift;
	
	# load the package
	my %map = (
		'clustalw'  => 'Clustalw',
		'kalign'    => 'Kalign',
		'mafft'     => 'MAFFT',
		'muscle'    => 'Muscle',
		'probalign' => 'Probalign',
		'probcons'  => 'Probcons',
		'tcoffee'   => 'TCoffee',
		'amap'      => 'Amap',
	);
	my $config  = $self->config;
	my $aligner = lc $config->MSA_TOOL || 'muscle';
	my $package = 'Bio::Tools::Run::Alignment::' . $map{$aligner};
	eval "require $package";
	$self->logger->error("couldn't load $package: $@") if $@;
	
	# instantiate
	my $instance = $package->new;
	
	# run post-processing
	my %pp = (
		'muscle'    => \&_muscle,
		'mafft'     => \&_quiet,
		'clustalw'  => \&_quiet,
		'kalign'    => \&_quiet,
		'probalign' => \&_verbose,
		'probcons'  => \&_verbose,
		'tcoffee'   => \&_quiet,
		'amap'      => \&_verbose,
	);
	$pp{$aligner}->($instance) if $pp{$aligner};
	return $instance;
}

=item profile_align_files

Aligns two sets of already aligned sequences (two file names) against each other, returns
a string representation of the alignment

=cut

sub profile_align_files {
    my ($self,$file1,$file2) = @_;   
    my $result = `muscle -profile -in1 $file1 -in2 $file2 -quiet`;
    return $result;
}

sub _temp_fasta {
    my $array = shift;
    my ( $fh, $filename ) = tempfile();
    for my $seq ( @{ $array } ) {
	print $fh '>', $seq->get_name, "\n", $seq->get_char, "\n";
    }
    close $fh;
    return $filename;
}

=back

=cut

1;