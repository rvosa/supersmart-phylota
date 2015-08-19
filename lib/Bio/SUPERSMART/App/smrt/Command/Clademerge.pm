package Bio::SUPERSMART::App::smrt::Command::Clademerge;

use strict;
use warnings;

use List::MoreUtils 'uniq';

use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse_matrix';
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Service::ParallelService;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: merges sets of alignments into input files for clade inference

=head1 NAME

Clademerge.pm - For each decomposed clade, merges the set of alignments assembled for this clade into an input file for tree inference.

=head1 SYNOPSYS

smrt clademerge [-h ] [-v ] [-w <dir>] [-o <format>]

=head1 DESCRIPTION

Given a working directory, traverses it looking for subdirectories of the pattern
clade*. Perusing each of these, it merges the *.fa (FASTA) files in it and
produces a single output file that can be analysed by the subcommand cladeinfer.

=cut


sub options {
    my ($self, $opt, $args) = @_;       
    return (
        [ "outformat|o=s", "output format for merged clade files (phylip or nexml), defaults to 'nexml'", { arg=>"format", default=> 'nexml'} ],
        [ "enrich|e", "enrich the selected markers with additional haplotypes", {} ],
    );  
}

sub validate {};

sub run {
    my ($self, $opt, $args) = @_;       
    
    my $workdir   = $self->workdir;
    my $outformat = $opt->outformat;
        
    # instantiate helper objects
    my $factory = Bio::Phylo::Factory->new;
    my $service = Bio::Phylo::PhyLoTA::Service::TreeService->new;
    my $mts     = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $sg      = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
    my $mt      = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $config  = Bio::Phylo::PhyLoTA::Config->new;
    my $log     = $self->logger;
    my $ns      = 'http://www.supersmart-project.org/terms#';

    # collect candidate dirs
    $log->info("Going to look for clade data in $workdir");
    my @dirs;
    opendir my $odh, $workdir or die $!;
    while( my $dir = readdir $odh ) {
        if ( $dir =~ /^clade\d+$/ and -d "${workdir}/${dir}" ) {
            push @dirs, $dir;
        }
    }

	# merge alignments for each clade 
	# This is done before the other steps, since SequenceGetter's merge_alignment
	# uses pmap and otherwise pmap calls would be nested.
	for my $cladedir ( @dirs ) {
		my $mergedfile = "${workdir}/${cladedir}/merged.txt";
		
		if ( -e $mergedfile  and -s $mergedfile) {
			$log->info("Alignments in clade $cladedir already merged");		
		}
		else {
			# collect seed gis from alignment file names
			my @files =  <"${workdir}/${cladedir}/*.fa">;
			my @gis =  grep { $_ ne '' } map {$1 if $_=~/\/([0-9]+)-clade[0-9]+\.fa/ } @files;
			
			# merge alignments
			$log->info("Merging alignments in clade $cladedir");		
			$sg->merge_alignments( $config->CLADE_MAX_DISTANCE, "${workdir}/${cladedir}", $mergedfile, @gis );
		}
	}
	
    # process in parallel
    my @result = grep { defined $_ and -e $_ } pmap {
        my $dir = $_;#( $dir ) = @_;
        
        # initialize the container objects
        my $project = $factory->create_project( '-namespaces' => { 'smrt' => $ns } );
        my $taxa    = $factory->create_taxa;
        $project->insert($taxa);        
        
        # start processing the directory
        $log->info("Going to enrich alignments in $workdir/$dir");
        my @matrices;

		# read list of merged alignment files
		my $mergedfile = "${workdir}/${dir}/merged.txt";
		$log->debug("Trying to open merged file $mergedfile");
		open my $fh, '<', $mergedfile or die $!;
		my @files;
		push @files, $_ while(<$fh>);
		chomp (@files);

		for my $file ( @files ) {
        
			# parse the file, enrich and degap it
			$log->info("Adding alignment $file");
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
			$log->info("Found more alignments in clade directory $dir than CLADE_MAX_MARKERS. Using the first " . $self->config->CLADE_MAX_MARKERS . " largest alignments.");
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
				$log->info("Removing taxon " . $taxname . " from $dir");
				$tax->delete($t);
				
				# remove rows from matrix containing the taxon
				for my $mat ( @matrices ) {
					for my $row ( @{$mat->get_entities} ) {
						if ( my $taxid = $row->get_name =~ /$taxname/ ) {
							$log->info("Removing row for taxon " . $taxname . " from matrix");
							$mat->delete($row);
						}
					}
				}
			}
		}
		
		# write table listing all marker accessions for taxa
		my @marker_table = $mts->get_marker_table( @{ $project->get_items(_MATRIX_) } );
		$mts->write_marker_table( "${workdir}/${dir}/${dir}-markers.tsv", \@marker_table );
		
        # write the merged nexml
        if ( lc $outformat eq 'nexml' ) {
		my $outfile = "${workdir}/${dir}/${dir}.xml";
		$log->info("Going to write file $outfile");
		open my $outfh, '>', $outfile or die $!;
		print $outfh $project->to_xml( '-compact' => 1 );
		return $outfile;
        }
        
        # write supermatrix phylip
        elsif ( lc $outformat eq 'phylip' ) {
		my $outfile = "${workdir}/${dir}/${dir}.phy";
		my @matrices = @{ $project->get_items(_MATRIX_) };
		my ($taxa) = @{ $project->get_items(_TAXA_) };
		$log->info("Going to write file $outfile");
		$service->make_phylip_from_matrix($taxa, $outfile, @matrices);
		return $outfile;
        }
    } @dirs;
    
    # report result
    $log->info("Wrote outfile $_") for @result;
    $log->info("DONE");
    return 1;
}

1;
