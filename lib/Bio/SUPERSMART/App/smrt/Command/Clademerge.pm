package Bio::SUPERSMART::App::smrt::Command::Clademerge;

use strict;
use warnings;

use List::MoreUtils 'uniq';

use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse_matrix';

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
        
    # instantiate helper objects
    my $mt      = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $config  = Bio::Phylo::PhyLoTA::Config->new;
    my $log     = $self->logger;


    # collect candidate dirs
    $log->info("Going to look for clade data in $workdir");
    my @dirs;
    opendir my $odh, $workdir or die $!;
    while( my $dir = readdir $odh ) {
        if ( $dir =~ /^clade\d+$/ and -d "${workdir}/${dir}" ) {
            push @dirs, $dir;
        }
    }

	# merge alignments for each clade, don't include this in the pmap below since
	# orthologize_cladedir also uses pmap, hence avoid nested pmap calls
	for my $cladedir ( @dirs ) {
		my $mergedfile = "${workdir}/${cladedir}/merged.txt";		
		$mt->orthologize_cladedir( 
			'dir'=>$cladedir, 
			'outfile'=>$mergedfile, 
			'maxdist'=>$config->CLADE_MAX_DISTANCE );
	}
	
    # enrich (if requested in argument) and write matrix to file
    my @result = grep { defined $_ and -e $_ } pmap {
		(my $clade) = @_;
		my $dir = "${workdir}/${clade}";
		$mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new("${dir}/merged.txt", $config->CLADE_MIN_COVERAGE);
		$mt->write_clade_matrix( 
			'markersfile' => "${dir}/${clade}-markers.tsv",
			'outfile' => $opt->outformat eq 'phylip' ? "${dir}/${clade}.phy" : "${dir}/${clade}.xml",
			'max_markers' => $config->CLADE_MAX_COVERAGE,
			'enrich' => $opt->enrich,
			'format' => $opt->outformat);
    } @dirs;
    
    # report result
    $log->info("Wrote outfile $_") for @result;
    $log->info("DONE");
    return 1;
}

1;
