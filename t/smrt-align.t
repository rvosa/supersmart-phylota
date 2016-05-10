#!/usr/bin/perl
use warnings;

BEGIN {
	our @plan = $ENV{'TEST_SMRT'} ? 'no_plan' : 'skip_all' => 'env var TEST_SMRT not set';
}
use Test::More @plan;

use App::Cmd::Tester;
use FindBin '$Bin';
use File::Temp qw(tempfile tempdir);
use Bio::SUPERSMART::App::smrt;

my $log;

BEGIN {
	use Bio::Phylo::Util::Logger ':levels';
	$log = Bio::Phylo::Util::Logger->new( 
		'-level' => WARN, 
		'-class' => [
			'main',
			'Bio::SUPERSMART::Service::ParallelService',
		]
	);
}

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

my $workdir = tempdir( CLEANUP => 1, 'DIR' => $ENV{HOME} );
my $outfile = $workdir . "/aligned.txt";
my $infile = "$Bin/testdata/species-fishes.tsv";

my $result = test_app( Bio::SUPERSMART::App::smrt=> [ "align",  "-w",   $workdir, "-i", $infile, "-o", $outfile ]);

_output_ok( $result, $outfile );
