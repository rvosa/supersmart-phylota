package Bio::SUPERSMART::App::smrt::Command::BBmerge;
use strict;
use warnings;
use Bio::Phylo::Factory;
use Bio::Phylo::Matrices::Datum;
use Bio::SUPERSMART::Config;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;
use Bio::SUPERSMART::Service::TreeService;
use Bio::SUPERSMART::App::SubCommand;
use List::MoreUtils qw(uniq);
use List::Util qw(max);
use Data::Dumper;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: creates supermatrix for genus-level backbone tree

=head1 NAME

BBmerge.pm - creates supermatrix for genus-level backbone tree

=head1 SYNOPSIS

smrt bbmerge [-h ] [-v ] [-w <dir>] -a <file> -t <file> [-o <file>] 

=head1 DESCRIPTION

Given an input file that lists alignment file locations (one on each line), traverses each 
genus in each alignment and picks the most divergent two species to represent 
their genus. The rationale is that these species will (likely) cross the root of their 
genus, so that the below genus-level tree can then be scaled to the same depth of that 
split and be grafted onto the tree without (too many) negative branch lengths.

The way in which the overall two divergent species within the genus are selected is as 
follows:

* for each alignment, within each genus, make all pairwise comparisons and sort
  the pairs by decreasing sequence divergence.

* pick the most distal pair and weight it in proportion to the number of pairs, 
  within that genus for that alignment, minus one. This means that singleton pairs are 
  discarded, and those from bigger samples are assumed to more accurately indicate which 
  taxa actually cross the root.

* after having processed all alignments, pick the species pair that has the highest
  score. 

Subsequently, the optimal combination of markers needs to be selected to best cover the
exemplars. It is not optimal to just concatenate all alignments that cover any of the 
taxa - this can result in monstrous, sparse, supermatrices. Instead we give the user the
possibility of assembling a set of alignments such that all exemplar species are covered
by at least some minimal value, (though relatively frequently studied species would exceed
this). This is done as follows:

* for each exemplar species, collect all alignments that include it and sort this 
  collection in decreasing exemplar taxon coverage (i.e. the first alignment has the most 
  exemplar species in it, the last alignment the fewest).

* sort the exemplar species by increasing overall participation in the alignments
  (i.e. the first exemplar has been sequenced the fewest times, the last one the most).

* iterate over the sorted list of exemplars, and for each exemplar add their not-yet-seen, 
  sorted alignments to the stack, one by one. After adding each alignment, update the 
  coverage counts for all exemplar species that participate in that alignment. End the 
  iterations when all exemplars have crossed their threshold or have no more alignments 
  available.

Point of consideration: the node depths on the exemplar tree will be underestimates
relative to the genus-level tree (due to the node density effect), so it might be better
to give the node depths from the exemplar tree as priors and then see where the posteriors
come out. Disadvantage is that this is likely to lead to negative branch lengths in some
cases.

=cut

sub options {
    my ( $self, $opt, $args ) = @_;
    my $outfile_default      = "supermatrix.phy";
    my $outformat_default    = "phylip";
    my $markerstable_default = "markers-backbone.tsv";
    my $taxa_default         = "species.tsv";
    my $merged_default       = "merged.txt";	
    my $config       = Bio::SUPERSMART::Config->new;
	my $exemplars_default = $config->BACKBONE_EXEMPLARS;
	my @formats = qw (bl2seq clustalw emboss fasta maf mase mega meme msf nexus pfam phylip prodom psi selex stockholm mrbayes);
    return (
        [
            "alnfile|a=s",
			"list of file locations of merged alignments  as produced by 'smrt orthologize'",
            { arg => "file", default => $merged_default, galaxy_in => 1, galaxy_format => 'tabular', galaxy_type => "data" }
        ],
        [
            "taxafile|t=s",
            "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'",
            { arg => "file", default => $taxa_default, galaxy_in => 1, galaxy_format => 'tabular', galaxy_type => "data" }
        ],
        [
            "outfile|o=s",
            "name of the output file, defaults to '$outfile_default'",
            { default => $outfile_default, arg => "file", galaxy_out => 1, galaxy_format => $outformat_default, galaxy_type => "data", galaxy_label => 'supermatrix' }
        ],
        [
            "format|f=s",
			"format of supermatrix, defaults to '$outformat_default'; possible formats: " . join (', ', @formats),
            { default => $outformat_default, galaxy_in => 1, galaxy_type => "select", galaxy_options => \@formats, galaxy_value => $outformat_default }
        ],
        [
            "include_taxa|i=s",
			"one or multiple names of taxa present in 'taxafile' (e.g. species or genus names, separated by commata) whose representative species will be included in the output dataset, regardless of marker coverage and sequence divergence",
            { galaxy_in => 1, galaxy_type => "text" }
        ],
        [
            "markersfile|m=s",
			"name for summary table with included accessions, defaults to $markerstable_default",
            { default => $markerstable_default, arg => "file", galaxy_out => 1, galaxy_format => 'tabular', galaxy_type => "data" }
        ],
        [
            "exemplars|e=s",
			"number of exemplar species per genus, defaults to $exemplars_default, set to -1 to include all species",
            { default => $exemplars_default, galaxy_in => 1, galaxy_type => "integer", galaxy_value => $exemplars_default }
        ],
        [
		 "names|n",
		 "print taxon names instead of identifiers in supermatrix",
		 { default => 0, galaxy_in => 1, galaxy_type => "boolean" }
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
    my $include_taxa = $opt->include_taxa;

    # instantiate helper objects
    my $config = Bio::SUPERSMART::Config->new;
    my $mt  = Bio::SUPERSMART::Domain::MarkersAndTaxa->new($opt->alnfile, $config->BACKBONE_MIN_COVERAGE);
    my $log = $self->logger;
        
    # Pick the exemplar taxa:
    # we first further narrow down the list of possible exemplars by the following 
    # criterion: A taxon must share at least one marker with a taxon in its own genus.
    # If after that we are still left with more than two candidates in a genus,
    # we pick the taxa for which the sequences have the highest divergence.
    my @exemplars = $mt->pick_exemplars( $taxafile, $include_taxa, $opt->exemplars );
    $log->info( "Identified " . scalar(@exemplars) . " exemplars" );

    # optimize the order in which taxa and alignments are added to the supermatrix.
    # this means starting at the least-sequenced taxon and the most speciose alignment
    my ($sorted_exemplars,$sorted_alignments) = $mt->optimize_packing_order(@exemplars);        
    $log->info( "Using " . scalar(@$sorted_alignments) . " alignments for supermatrix" );
    $log->info( "Number of exemplars : " . scalar(@$sorted_exemplars) );

    # write alignmnets to supermatrix file
    $mt->write_supermatrix( 
    	'alignments'  => $sorted_alignments, 
    	'exemplars'   => $sorted_exemplars,
        'outfile'     => $self->outfile, 
        'format'      => $opt->format, 
        'markersfile' => $opt->markersfile, 
		'taxon_names'       => $opt->names,
    );
    $log->info("DONE, results written to $outfile");
    return 1;
}

1;
