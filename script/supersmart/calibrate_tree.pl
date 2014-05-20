#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::CalibrationService;
use Bio::Phylo::PhyLoTA::Service::TreeService;

=head1 NAME

calibrate_tree.pl - time calibrates a tree using treePL

=head1 SYNOPSYS

 $ calibrate_tree.pl --fossiltable=<file> --readtree=<file> \
 	--numsites=<number of sites> --smooth=<PL smoothing value> > --outfile=<file>
 
=head1 DESCRIPTION

Given a spreadsheet of fossils and an input tree, runs the treePL program to calibrate the tree using penalized likelihood. 
The treePL algorithm requires a smoothing parameter and number of sites, which must be provided as
arguments. Also the number of sites in the alignment has to be provided. Writes the calibrated tree in newick format to the 
file specified by the 'outfile' argument. Also, a second newick file is produced with the taxon ids back to the taxon named. 
This file is named labelled_<outfile> by default.

=cut

# process command line arguments
my $fossiltable;  # tab-separated spreadsheet, as in data/FOSSILS/
my $readtree;     # the maximum likelihood tree to smooth
my $numsites;     # number of sites in the alignment that resulted in $readtree
my $outfile;    # the output tree to be written by treePL
my $smooth = 100; # PL smoothing value
my $verbosity = WARN; # logging level
GetOptions(
	'fossiltable=s' => \$fossiltable,
	'readtree=s'    => \$readtree,
	'verbose+'      => \$verbosity,
	'numsites=i'    => \$numsites,
	'outfile=s'   => \$outfile,
	'smooth=i'      => \$smooth,
    );

# set up helper objects
Bio::Phylo::Util::Logger->VERBOSE( '-level' => $verbosity, '-class' => [qw (main Bio::Phylo::PhyLoTA::Service::CalibrationService)] );
my $config = Bio::Phylo::PhyLoTA::Config->new;
my $cs = Bio::Phylo::PhyLoTA::Service::CalibrationService->new;
my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;

# read the ML tree from newick (one tree in this file)
INFO "reading tree from file $readtree";
my $tree = parse_tree( 
        '-format'     => 'newick', 
        '-file'       => $readtree, 
        '-as_project' => 1 
    );

my $table;
if ( $fossiltable ) {
        # read the spreadsheet with calibration points
        INFO "reading fossils from file $fossiltable";

        my @fossils = $cs->read_fossil_table($fossiltable);

        # make calibration table from fossils
        INFO "going to make calibration table";
        my @points = map { $cs->find_calibration_point($_) } @fossils;
        $table = $cs->create_calibration_table( $tree, @points );

}

# calibrate the tree
my $chronogram = $cs->calibrate_tree (('-numsites' => $numsites, '-calibration_table' => $table, '-tree' => $tree));

# save calibrated tree to file
open my $fh, '>', $outfile or die $!;
print $fh $chronogram->to_newick."\n";
close $fh;

# remap calibrated tree and save to file
my $labelled_chronogram = $ts->remap_to_name($chronogram);
my ( $volume, $directories, $file ) = File::Spec->splitpath( $outfile );
my $filename_lab = $volume .  $directories . 'labelled_' . $file;
open my $fh_lab, '>', $filename_lab or die $!;
print $fh_lab $labelled_chronogram->to_newick."\n";
close $fh_lab;
