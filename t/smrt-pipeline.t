#!/usr/bin/perl
use warnings;
BEGIN {
	our @plan = $ENV{'TEST_SMRT'} ? 'no_plan' : 'skip_all' => 'env var TEST_SMRT not set';
}
use Test::More @plan;
use App::Cmd::Tester;
use FindBin '$Bin';
use File::Temp qw(tempdir);

use Bio::Phylo::Util::Logger ':levels';
use Bio::SUPERSMART::App::smrt;

BEGIN { 
	use_ok('Bio::SUPERSMART::App::smrt::Command::Taxize'); 
	use_ok('Bio::SUPERSMART::App::smrt::Command::Align');
	use_ok('Bio::SUPERSMART::App::smrt::Command::Orthologize');
	use_ok('Bio::SUPERSMART::App::smrt::Command::BBreroot');
	use_ok('Bio::SUPERSMART::App::smrt::Command::BBmerge');
	use_ok('Bio::SUPERSMART::App::smrt::Command::BBinfer');
	use_ok('Bio::SUPERSMART::App::smrt::Command::BBcalibrate');
	use_ok('Bio::SUPERSMART::App::smrt::Command::BBdecompose');
	use_ok('Bio::SUPERSMART::App::smrt::Command::Clademerge');
	use_ok('Bio::SUPERSMART::App::smrt::Command::Cladeinfer');
	use_ok('Bio::SUPERSMART::App::smrt::Command::Cladegraft');
};

my $namesfile = "$Bin/../examples/fishes/names.txt";
my $fossilfile = "$Bin/testdata/marine_calibration.txt";

my $log = Bio::Phylo::Util::Logger->new( '-level' => INFO, '-class' => 'main' );

my $workdir = tempdir( CLEANUP => 0, 'DIR' => $ENV{HOME} );
my $taxafile = "taxa.tsv";
my $alnfile = "aln.txt";
my $mergefile = "merged.txt";
my $smfile = "supermatrix.phy";
my $bbfile = "backbone.dnd";
my $bbrefile = "bbrerooted.dnd";
my $chrofile = "chronogram.dnd";
my $summfile = "summary.txt";
my $finalfile = "final.dnd";

# test if run gave no errors and output file is good
sub _output_ok {
	my ( $cmd, $result, $outfile ) = @_;
	if ( $result->error ){
		print "ERROR :" .  $result->error . "\n";
	}
	ok ( $result->output =~ /DONE/, "command $cmd prints 'DONE'" );
	ok ( ! $result->error, "command $cmd : no error in result" );
	is ( $result->run_rv, 1, "command $cmd : return value is 1" );
	if ( $outfile ) {
		$outfile = $workdir . "/" . $outfile;
		ok ( -e $outfile, "outfile $outfile exists" );
		ok ( -s $outfile, "outfile $outfile not empty" );
	}
}

my $result = test_app( Bio::SUPERSMART::App::smrt=> [ "taxize",  "-i", $namesfile, "-o", $taxafile, "-e", "Species", "-w", $workdir, "-l", "taxize.log" ]);
_output_ok ( "taxize", $result, $taxafile );

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "align",  "-i", $taxafile, "-o", $alnfile , "-w", $workdir , "-l", "align.log" ]);
_output_ok ( "align", $result, $alnfile );

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "orthologize",  "-i",   $alnfile, "-o", $mergefile , "-w", $workdir , "-l", "orthologize.log" ]);
_output_ok ( "orthologize", $result, $mergefile );

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "bbmerge",  "-a",   $mergefile, "-o", $smfile, "-t", $taxafile , "-w", $workdir , "-l", "bbmerge.log" ]);
_output_ok ( "bbmerge", $result, $smfile );

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "bbinfer",  "-s",   $smfile, "-o", $bbfile, "-t", $taxafile, "-w", $workdir , "-i", "examl", "-l", "bbinfer.log" ]);
_output_ok ( "bbinfer", $result, $bbfile );

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "bbreroot",  "-b",  $bbfile, "-o", $bbrefile, "-t", $taxafile , "-w", $workdir , "-l", "bbreroot.log" ]);
_output_ok ( "bbreroot", $result, $bbrefile );

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "bbcalibrate",  "-s",   $smfile, "-f", $fossilfile, "-t", $bbrefile, "-o", $chrofile , "-w", $workdir , "-l", "bbcalibrate.log" ]);
_output_ok ( "bbcalibrate", $result, $chrofile );

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "bbdecompose",  "-b",  $chrofile, "-o", $summfile, "-t", $taxafile, "-a", $alnfile , "-w", $workdir , "-l", "bbdecompose.log" ]);
_output_ok ( "bbdecompose", $result, $summfile );

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "clademerge" , "-w", $workdir , "-l", "clademerge.log" ]);
_output_ok ( "clademerge", $result );

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "cladeinfer" , "-w", $workdir , "-l", "cladeinfer.log" ]);
_output_ok ( "cladeinfer", $result );

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "cladegraft",  "-b",  $chrofile, "-o", $finalfile , "-w", $workdir , "-l", "cladegraft.log" ]);
_output_ok ( "cladegraft", $result, $finalfile );



