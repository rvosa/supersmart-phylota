package Bio::SUPERSMART::App::smrt::Command::BBcalibrate;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Service::CalibrationService;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::ParallelService;

use Bio::Phylo::IO qw(parse_tree);

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);


# ABSTRACT: calibrate backbone tree using TreePL

=head1 NAME

BBcalibrate.pm - inference of genus-level backbone tree

=head1 SYNOPSIS

smrt bbcalibrate [-h ] [-v ] [-w <dir>] -t <file> -s <file> -f <file> [-o <file>] 

=head1 DESCRIPTION

Given a rooted molecular backbone phylogeny and a table of fossils, 
performs tree calibration using treePL. Produces an ultrametric tree 
with branch lengths proportional to evolutionary time (i.e. a "chronogram").

Runs the treePL program to calibrate the tree using penalized likelihood. 
The treePL algorithm requires a smoothing parameter and number of sites, 
which must be provided as arguments. Also the supermatrix which has been 
used for tree inference has to be provided because the number of alignment 
sites has to be given to TreePL. Writes the calibrated tree in newick 
format to the file specified by the 'outfile' argument. 

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "chronogram.dnd";
	my $tree_default    = "backbone-rerooted.dnd";
	my $matrix_default  = "supermatrix.phy";
	return (
		["fossiltable|f=s", "tsv (tab-separated value) file containing fossil table with at least 5 columns (id, name, crown/stem, taxon, age)", { arg => "file", mandatory => 1}],
		["tree|t=s", "backbone tree to calibrate as produced by 'smrt bbreroot'", { arg => "file", default => $tree_default}],
		["supermatrix|s=s", "matrix of concatenated multiple sequece alignments which was used to generate the tree", { arg => "file", default => $matrix_default}],	
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
	
	# prepare reusable variables
	$logger->info( "Reading fossils from file $fossiltable" );
	my @fossils = $cs->read_fossil_table($fossiltable);
	my @points = map { $cs->find_calibration_point($_) } @fossils;
	my $numsites = $mt->get_supermatrix_numsites($supermatrix);
	
	# read and write trees one by one
	open my $in, '<', $treefile or die $!;
	chomp(my @backbone_trees = <$in>);
	close $in;

	# mapping tables for faster id and taxon name lookup
	my %ti_to_name;
	my %name_to_ti;
	
	my @calibrated_trees = pmap {
		my $newick = $_;
		
		$logger->debug("Attempting to calibrate tree: $newick");
		my $tree = parse_tree( 
			'-format' => 'newick', 
			'-string' => $newick, 
			);
		$tree->resolve;

		# create id mapping table
		if ( ! scalar(%ti_to_name) ) {
			%ti_to_name = $ts->make_mapping_table($tree);
			%name_to_ti = reverse(%ti_to_name);
		}

		# map identifiers
		$tree = $ts->remap($tree, %name_to_ti);
		my $treestr = $tree->to_newick;
		
		# make calibration table from fossils
		$logger->info( "Going to make calibration table" );
		my $table = $cs->create_calibration_table( $tree, @points );
		if (! $table ) {
			$logger->warn("Could not create calibration table");
			return;
		}
		
		# refresh tree, it can happen that it gets damaged. This is a workaround
		$tree = parse_tree( 
			'-format' => 'newick', 
			'-string' => $treestr, 
			);

		# calibrate the tree
		my $nthreads = int($config->NODES/scalar(@backbone_trees)) || 1;
		my $chronogram = $cs->calibrate_tree (
			'-numsites'          => $numsites, 
			'-calibration_table' => $table, 
			'-tree'              => $tree,
			'-nthreads'          => $nthreads
			);
		
		# translate from taxon id's to names
		my $labelled_chronogram = $ts->remap($chronogram, %ti_to_name);

		$logger->info("Finished calibrating single tree");
		#	return $labelled_chronogram->to_newick;
		return $labelled_chronogram->to_newick;
	} @backbone_trees;
	
	$logger->warn("Number of calibrated trees different than number of input trees") if scalar(@backbone_trees) != scalar(@calibrated_trees);
	
	# write output file
	if ( my $cnt = scalar @calibrated_trees ) {
		$logger->info("Writing $cnt calibrated trees to $outfile");
		open my $out, '>', $outfile or die $!;			
		print $out $_  . "\n" for @calibrated_trees;
		close $out;
	}  
	else {
		$logger->warn("No calibrated trees to write!");
	}
	
	$logger->info("DONE, results written to $outfile");
	
	return 1;		
}

1;
