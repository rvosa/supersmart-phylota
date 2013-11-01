#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::FossilDataGetter;
use Bio::Phylo::PhyLoTA::Service::CalibrationTableCreator;

# process command line arguments
my $fossiltable;  # tab-separated spreadsheet, as in data/FOSSILS/
my $readtree;     # the maximum likelihood tree to smooth
my $numsites;     # number of sites in the alignment that resulted in $readtree
my $writetree;    # the output tree to be written by treePL
my $smooth = 100; # PL smoothing value
my $verbosity = WARN; # logging level
GetOptions(
	'fossiltable=s' => \$fossiltable,
	'readtree=s'    => \$readtree,
	'verbose+'      => \$verbosity,
	'numsites=i'    => \$numsites,
	'writetree=s'   => \$writetree,
	'smooth=i'      => \$smooth,
);

# configure logging
Bio::Phylo::Util::Logger->VERBOSE( 
	'-level' => $verbosity, 
	'-class' => 'main',
);

# read the ML tree from newick (one tree in this file)
INFO "reading tree from file $readtree";
my $tree = parse_tree( 
	'-format'     => 'newick', 
	'-file'       => $readtree, 
	'-as_project' => 1 
);

# read the spreadsheet with calibration points
INFO "reading fossils from file $fossiltable";
my $fdg = Bio::Phylo::PhyLoTA::Service::FossilDataGetter->new;
my @fossils = $fdg->read_fossil_table($fossiltable);

# make calibration table from fossils
INFO "going to make calibration table";
my $ctc = Bio::Phylo::PhyLoTA::Service::CalibrationTableCreator->new;
my @points = map { $ctc->find_calibration_point($_) } @fossils;
my $table = $ctc->create_calibration_table( $tree, @points );

# write the config file header
INFO "printing treePL config file header";
print <<"HEADER";
treefile = $readtree
smooth = $smooth
numsites = $numsites
outfile = $writetree
HEADER

# write the mrca statements
my $counter;
for my $row ( $table->get_rows ) {
	my @taxa = $row->taxa;
	print "mrca = xx" . ++$counter . " @taxa\n";
	print "max = xx$counter " . $row->max_age . "\n";
	print "min = xx$counter " . $row->min_age . "\n";
	INFO "printed mrca $counter for @taxa";
	INFO "minimum age is " . $row->min_age;
	INFO "maximum age is " . $row->max_age;
}

# we do a "thorough" analysis, see the treePL wiki on github
print "thorough\n";

