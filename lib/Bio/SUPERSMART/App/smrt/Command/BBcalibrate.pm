package Bio::SUPERSMART::App::smrt::Command::BBcalibrate;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Service::CalibrationService;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Config;

use Bio::Phylo::IO qw(parse_tree);

use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);


# ABSTRACT: calibrate backbone tree using TreePL

=head1 NAME

BBcalibrate.pm - inference of genus-level backbone tree

=head1 SYNOPSIS

smrt bbcalibrate [-h ] [-v ] [-w <dir>] -t <file> -s <file> -f <file> [-o <file>] 

=head1 DESCRIPTION

Given a rooted molecular backbone phylogeny and a table of fossils, performs tree calibration using treePL. 
Produces an ultrametric tree with branch lengths proportional to evolutionary time (i.e. a "chronogram").

Runs the treePL program to calibrate the tree using penalized likelihood. 
The treePL algorithm requires a smoothing parameter and number of sites, which must be provided as
arguments. Also the supermatrix which has been used for tree inference has to be provided because the number of alignment sites has
to be given to TreePL. Writes the calibrated tree in newick format to the 
file specified by the 'outfile' argument. 

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "chronogram.dnd";
	return (
		["tree|t=s", "backbone tree to calibrate as produced by 'smrt bbinfer' or 'smrt bbreroot'", { arg => "file", mandatory => 1}],
		["supermatrix|s=s", "matrix of concatenated multiple sequece alignments which was used to generate the tree", { arg => "file", mandatory => 1}],	
		["fossiltable|f=s", "tsv (tab-separated value) file containing fossil table with at least 5 columns (id, name, crown/stem, taxon, age)", { arg => "file", mandatory => 1}],	
		["outfile|o=s", "name of the output tree file (in newick format), defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],			

	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	#  Abort if input files not present or empty  
	my @files = ( $opt->supermatrix, $opt->tree, $opt->fossiltable );
	foreach my $file ( @files ){
		$self->usage_error("need tree, supermatrix and fossiltable arguments") if not $file;
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
}

sub run {
	my ($self, $opt, $args) = @_;		
	
	# collect command-line arguments
	my $treefile = $opt->tree;
	my $supermatrix = $opt->supermatrix;
	my $fossiltable = $opt->fossiltable;
	my $outfile = $self->outfile;
	my $logger = $self->logger;
	 
	# instantiate helper objects
	my $config = Bio::Phylo::PhyLoTA::Config->new;
	my $cs = Bio::Phylo::PhyLoTA::Service::CalibrationService->new;
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	
	# read the ML tree from newick (one tree in this file)
	$logger->info( "reading tree from file $treefile" );
	my $tree = parse_tree( 
        '-format'     => 'newick', 
        '-file'       => $treefile, 
        '-as_project' => 1 
    );
	$tree->resolve;
	$tree = $ts->remap_to_ti($tree);
	
    # read the table with calibration points
    $logger->info( "reading fossils from file $fossiltable" );

    my @fossils = $cs->read_fossil_table($fossiltable);

    # make calibration table from fossils
    $logger->info( "going to make calibration table" );
    my @points = map { $cs->find_calibration_point($_) } @fossils;
    my $table = $cs->create_calibration_table( $tree, @points );

	my $numsites = $mt->get_supermatrix_numsites($supermatrix);

	# calibrate the tree
	my $chronogram = $cs->calibrate_tree (('-numsites' => $numsites, '-calibration_table' => $table, '-tree' => $tree));
	
	# translate from taxon id's to names and save to file
	my $labelled_chronogram = $ts->remap_to_name($chronogram);
	open my $outfh, '>', $outfile or die $!;
	print $outfh $labelled_chronogram->to_newick;
	close $outfh;
	
	$logger->info("DONE, results written to $outfile");
	
	return 1;		
}

1;