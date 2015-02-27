#!/usr/bin/perl
use warnings;

use Test::More 'no_plan';
use App::Cmd::Tester;
use FindBin '$Bin';
use File::Temp qw(tempfile tempdir);

use Bio::SUPERSMART::App::smrt;

BEGIN { use_ok('Bio::SUPERSMART::App::smrt::Command::BBinfer') };

# test if run gave no errors and output file is good
sub _output_ok {
	my ( $result, $outfile ) = @_;
	ok ( $result->output =~ /DONE/, "command prints 'DONE'" );
	ok ( ! $result->error, "no error in result" );
	is ( $result->run_rv, 1, "return value is 1" );
	ok ( -e $outfile, "outfile exists" );
	ok ( -s $outfile, "outfile not empty" );
}

my $matrix = "$Bin/testdata/supermatrix-fishes.phy";
my $classtree = "$Bin/testdata/classificationtree-fishes.dnd";

my ($fh, $outfile) = tempfile( 'CLEANUP' => 1 );

my $result = test_app( Bio::SUPERSMART::App::smrt=> [ "bbinfer",  "-s",   $matrix, "-t", $classtree, "-o", $outfile ]);
_output_ok ( $result, $outfile );

# test inference with raxml and bootstrapping
$result = test_app( Bio::SUPERSMART::App::smrt=> [ "bbinfer",  "-s",   $matrix, "-t", $classtree, "-i", "raxml", "-o", $outfile, "-b"]);
_output_ok ( $result, $outfile );

# test inference with exabayes
$result = test_app( Bio::SUPERSMART::App::smrt=> [ "bbinfer",  "-s",   $matrix, "-i", "exabayes", "-o", $outfile, "-b"]);
_output_ok ( $result, $outfile );
