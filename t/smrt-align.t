#!/usr/bin/perl
use warnings;

use Test::More 'no_plan';
use App::Cmd::Tester;
use FindBin '$Bin';
use File::Temp qw(tempfile tempdir);

use Bio::Phylo::Util::Logger ':levels';
use Bio::SUPERSMART::App::smrt;

my $log = Bio::Phylo::Util::Logger->new( '-level' => INFO, '-class' => 'main' );

BEGIN { use_ok('Bio::SUPERSMART::App::smrt::Command::Align') };

# test if run gave no errors and output file is good
sub _output_ok {
	my ( $result, $outfile ) = @_;
	ok ( $result->output =~ /DONE/, "command prints 'DONE'" );
	if ( $result->error ){
		$log->info($result->error);
	}
	ok ( ! $result->error, "no error in result" );
	is ( $result->run_rv, 1, "return value is 1" );
	ok ( -e $outfile, "outfile exists" );
	ok ( -s $outfile, "outfile not empty" );
}

my $workdir = tempdir( CLEANUP => 0 );
my $outfile = $workdir . "/aligned.txt";
my $infile = "$Bin/testdata/species-fishes.tsv";

$log->info("Workdir : $workdir");

my $result = test_app( Bio::SUPERSMART::App::smrt=> [ "align",  "-w",   $workdir, "-i", $infile, "-o", $outfile ]);

_output_ok( $result, $outfile );