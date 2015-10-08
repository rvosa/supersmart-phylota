package Bio::SUPERSMART::App::smrtutils::Command::Supermatrix;
use strict;
use warnings;
use Bio::Phylo::Factory;
use Bio::Phylo::Matrices::Datum;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Config;
use Bio::SUPERSMART::App::SubCommand;
use List::MoreUtils qw(uniq);
use List::Util qw(max);
use Data::Dumper;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: 

=head1 NAME



=head1 SYNOPSIS



=head1 DESCRIPTION


=cut

sub options {
    my ( $self, $opt, $args ) = @_;
    my $outfile_default      = "supermatrix-all.phy";
    my $outformat_default    = "phylip";
    my $markerstable_default = "markers-supermatrix.tsv";
    my $taxa_default         = "species.tsv";
    my $merged_default       = "merged.txt";	
    my $config       = Bio::Phylo::PhyLoTA::Config->new;
	my $exemplars_default = $config->BACKBONE_EXEMPLARS;
    return (
        [
            "alnfile|a=s",
			"list of file locations of merged alignments  as produced by 'smrt orthologize'",
            { arg => "file", default => $merged_default }
        ],
        [
            "taxafile|t=s",
            "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'",
            { arg => "file", default => $taxa_default }
        ],
        [
            "outfile|o=s",
            "name of the output file, defaults to '$outfile_default'",
            { default => $outfile_default, arg => "file" }
        ],
        [
            "format|f=s",
			"format of supermatrix, defaults to '$outformat_default'; possible formats: phylip, nexml",
            { default => $outformat_default }
        ],
        [
            "markersfile|m=s",
			"name for summary table with included accessions, defaults to $markerstable_default",
            { default => $markerstable_default, arg => "file" }
        ],
        [
		    "enrich|r",
			"enrich the selected markers with additional haplotypes",
		 {}
        ],
    );
}

sub validate {
    my ( $self, $opt, $args ) = @_;

    # If alignment or taxa file is absent or empty, abort
    my @files = ( $opt->alnfile, $opt->taxafile );
    for my $file (@files) {
        $self->usage_error("need alnfile and taxafile arguments") if not $file;
        $self->usage_error("file $file does not exist") unless -e $file;
        $self->usage_error("file $file is empty")       unless -s $file;
    }
}

sub run {
    my ( $self, $opt, $args ) = @_;

    # collect command-line arguments
    my $taxafile     = $opt->taxafile;
    my $outfile      = $self->outfile;
	
    # instantiate helper objects
    my $mt  = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new($opt->alnfile);
    my $log = $self->logger;
	my $config  = Bio::Phylo::PhyLoTA::Config->new;

    # Pick all species with sufficient data
    my @species = $mt->pick_exemplars( $taxafile, 0, -1 );
	$log->info( scalar(@species) . ' species have sufficient data for supermatrix');

	my $filename = $mt->write_clade_matrix( 
		
		'outfile'     => $outfile,
		'markersfile' => $opt->markersfile,
		'min_markers' => $config->CLADE_TAXON_MIN_MARKERS,
		'max_markers' => $config->CLADE_MAX_MARKERS,
		'enrich' => $opt->enrich,
		'format' => $opt->format);

    $log->info("DONE, results written to $outfile");
    return 1;
}

1;
