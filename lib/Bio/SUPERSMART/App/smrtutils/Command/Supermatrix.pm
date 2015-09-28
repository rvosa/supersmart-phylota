package Bio::SUPERSMART::App::smrtutils::Command::Supermatrix;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::Factory;
use Bio::Phylo::PhyLoTA::Config;

use Bio::Phylo::IO qw(parse);
use Bio::Phylo::Util::CONSTANT ':objecttypes';

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: Creates supermatrix for all taxa

=head1 NAME

Supermatrix.pm.pm - Compiles one supermatrix for all taxa containing all data

=head1 SYNOPSYS

smrt-utils 

=head1 DESCRIPTION

=cut

sub options {    
	my ($self, $opt, $args) = @_;
	my $outfile_default = "supermatrix-alltaxa.phy";
	my $alnfile_default = "aligned.txt";
	my $taxafile_default = "species.tsv";
	return (
		['outfile|o=s', "output file, defaults to $outfile_default", { arg => 'file', default => $outfile_default }],    	    
		['alnfile|a=s', "alignment file, defaults to $alnfile_default", { arg => 'file', default => $alnfile_default }],    	    		
		['taxafile|t=s', "taxa file, defaults to $taxafile_default", { arg => 'file', default => $taxafile_default }],    	    		
		);	
}

sub validate {
    my ($self, $opt, $args) = @_;
}

sub run {
	my ($self, $opt, $args) = @_;    
	my $logger = $self->logger;

	my $outfile = $opt->outfile;
	my $alnfile = $opt->alnfile;
	my $taxafile = $opt->taxafile;
	my $workdir = $self->workdir;

    my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
    my $factory = Bio::Phylo::Factory->new;
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $sg  = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
    my $config = Bio::Phylo::PhyLoTA::Config->new;
    my $service = Bio::Phylo::PhyLoTA::Service::TreeService->new;
    my $ns      = 'http://www.supersmart-project.org/terms#';

	# read taxa table and get taxon IDs
	my @taxa = $mt->parse_taxa_file($taxafile);
	my @taxids = map {$_->{'species'}} @taxa;

	# read alignment list 
	my @alignments;
    open my $fh, '<', $alnfile or die $!;
    while(<$fh>) {
        chomp;
        push @alignments, $_ if /\S/ and -e $_;
    }
	
	# filter for density and distance 
	my @filtered_alns = $mts->filter_clade_alignments( '-ingroup' => \@taxids, '-alnfiles' => \@alignments);
	$logger->info("Filtered out " . scalar(@filtered_alns) . " alignments");

	# determine the largest subset of species connected by markers
	my %adj = $mts->generate_marker_adjacency_matrix(\@filtered_alns, \@taxids);
	my @subsets = @{$mts->get_connected_taxa_subsets(\%adj)};
	
	# sort subsets by size, decreasing ans select largest one
	@subsets = sort { scalar(@$b) <=> scalar(@$a)  } @subsets;	
	my %aln_taxa = map {$_=>1} @{$subsets[0]};

	# write filtered alignents in temporary directory
	my $tmpdir = "${workdir}/tmp-supermatrix";
	mkdir( $tmpdir );
	my @to_merge;
	for my $al ( @filtered_alns ) {		
		my $seed_gi;
		my $seqstr;
		for my $k ( keys %$al ) {
			$seed_gi = $1 if $k =~ /seed_gi\|([0-9]+)/;
			my $ti = $1 if $k =~ /taxon\|([0-9]+)/;			
			$seqstr .= "> $k\n" . $al->{$k} . "\n" if $aln_taxa{$ti};
		}
		my $r = int(rand(1000)); # there can be duplicate seed gis, add random number to filename
		my $filename = "${tmpdir}/${seed_gi}-tmpaln-${r}.fa";
		push @to_merge, $filename;
		open my $fh, '>', $filename or die $!;
		print $fh $seqstr;
		close $fh;
	}
	
    # merge alignments
	my $mergedfile = "${tmpdir}/merged.txt";
	
	
	# Below, copied code from Clademerge.pm. TODO: Declutter clademerge, put this code
	# into service class
	my @gis =  grep { $_ ne '' } map {$1 if $_=~/\/([0-9]+)-tmpaln-[0-9]+\.fa/ } @to_merge;
	$sg->merge_alignments( $config->CLADE_MAX_DISTANCE, $tmpdir, $mergedfile, @gis );
	
	# initialize the container objects
	my $project = $factory->create_project( '-namespaces' => { 'smrt' => $ns } );
	my $taxa    = $factory->create_taxa;
	$project->insert($taxa);        
	
	# start processing the directory
	$logger->info("Going to enrich alignments in $tmpdir");
	my @matrices;
	
	# read list of merged alignment files			
	$logger->debug("Trying to open merged file $mergedfile");
	open my $ffh, '<', $mergedfile or die $!;
	my @files;
	push @files, $_ while(<$ffh>);
	chomp (@files);
	
	for my $file ( @files ) {
        
		# parse the file, enrich and degap it
		$logger->info("Adding alignment $file");
		my $matrix = $mt->parse_fasta_as_matrix(
			'-name' => $file,
			'-file' => $file,
			'-taxa' => $taxa,
			);
		$mts->enrich_matrix($matrix) if $opt->enrich;
		$matrix = $mts->degap_matrix($matrix);
		push @matrices, $matrix if $matrix;
		
	}
	return undef if not @matrices;
	
	# pick CLADE_MAX_MARKERS biggest alignments
	@matrices = sort { $b->get_ntax <=> $a->get_ntax } @matrices;
	if ( scalar(@matrices) > $self->config->CLADE_MAX_MARKERS ) {
		$logger->info("Found more alignments in clade directory $tmpdir than CLADE_MAX_MARKERS. Using the first " . $self->config->CLADE_MAX_MARKERS . " largest alignments.");
	}
	
	# keep track of taxon ids, the number of markers for each taxid,
	# because some markers might not be selected and taxa have to be removed
	my %markers_for_taxon;
	for my $i ( 0 .. $self->config->CLADE_MAX_MARKERS - 1 ) {
		if ( my $mat = $matrices[$i] ) {
			my @ids_for_mat;
			for ( @{ $mat->get_entities } ) {
					my $taxid = $_->get_name;
					$taxid =~ s/_.$//g;
					push @ids_for_mat, $taxid;
			}
			$markers_for_taxon{$_}++ for uniq (@ids_for_mat);
			$project->insert($mat);
			}
	}
	
	# remove a taxon from matrix if it has none or less markers than given in CLADE_TAXON_MIN_MARKERS
	# also remove all rows in the matrices where this taxon appears
	my ($tax) = @{ $project->get_items(_TAXA_) } ;
	for my $t ( @ {$tax->get_entities} ) {
		my $taxname = $t->get_name;
		my $marker_cnt = $markers_for_taxon{$taxname};
		if ( ! $marker_cnt || $marker_cnt < $self->config->CLADE_TAXON_MIN_MARKERS ) {
			$logger->info("Removing taxon " . $taxname . " from $tmpdir");
			$tax->delete($t);
			
			# remove rows from matrix containing the taxon
			for my $mat ( @matrices ) {
				for my $row ( @{$mat->get_entities} ) {
					if ( my $taxid = $row->get_name =~ /$taxname/ ) {
							$logger->info("Removing row for taxon " . $taxname . " from matrix");
							$mat->delete($row);
					}
				}
			}
		}
		}
	
	# write table listing all marker accessions for taxa
	#my @marker_table = $mts->get_marker_table( @{ $project->get_items(_MATRIX_) } );
	#$mts->write_marker_table( "${workdir}/${dir}/${dir}-markers.tsv", \@marker_table );
		
	# write supermatrix phylip
	
	#y $outfile = "${workdir}/${dir}/${dir}.phy";
	#my @matrices = @{ $project->get_items(_MATRIX_) };
	#my ($taxa) = @{ $project->get_items(_TAXA_) };
	#$logger->info("Going to write file $outfile");
	#$service->make_phylip_from_matrix($taxa, $outfile, @matrices);
		


}

1;
