package Bio::SUPERSMART::App::smrt::Command::Clademerge;

use strict;
use warnings;

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
    my $mt      = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
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
    
    # process in parallel
    my @result = grep { defined $_ and -e $_ } pmap {
        my ( $dir ) = @_;
        
        # initialize the container objects
        my $project = $factory->create_project( '-namespaces' => { 'smrt' => $ns } );
        my $taxa    = $factory->create_taxa;
        $project->insert($taxa);        
        
        # start processing the directory
        $log->info("Going to merge alignments in $workdir/$dir");
        my @matrices;
        opendir my $idh, "${workdir}/${dir}" or die $!;
        while( my $file = readdir $idh ) {
        
            # found an alignment
            if ( $file =~ /(.+?)\.fa$/ ) {
                my $id = $1;
                
                # parse the file, enrich and degap it
                $log->info("Adding alignment $file");
                my $matrix = $mt->parse_fasta_as_matrix(
                    '-name' => $id,
                    '-file' => "${workdir}/${dir}/${file}",
                    '-taxa' => $taxa,
                );
                $mts->enrich_matrix($matrix) if $opt->enrich;
                $matrix = $mts->degap_matrix($matrix);
                push @matrices, $matrix if $matrix;
            }
        }
        return undef if not @matrices;
        
        # pick CLADE_MAX_MARKERS biggest alignments
        @matrices = sort { $b->get_ntax <=> $a->get_ntax } @matrices;
		if ( scalar(@matrices) > $self->config->CLADE_MAX_MARKERS ) {
			$log->info("Found more alignments in clade directory $dir than CLADE_MAX_MARKERS. Using the first " . $self->config->CLADE_MAX_MARKERS . " alignments.");
		}
        for my $i ( 0 .. $self->config->CLADE_MAX_MARKERS - 1 ) {
            $project->insert($matrices[$i]) if $matrices[$i];
        }
        
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
            my ($taxa) = @{$project->get_items(_TAXA_)} ;
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
