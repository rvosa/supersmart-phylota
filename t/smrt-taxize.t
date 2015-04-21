#!/usr/bin/perl
use warnings;

use Test::More 'no_plan';
use App::Cmd::Tester;
use FindBin '$Bin';
use File::Temp qw(tempfile tempdir);

use Bio::SUPERSMART::App::smrt;

BEGIN { use_ok('Bio::SUPERSMART::App::smrt::Command::Taxize') };

# test if run gave no errors and output file is good
sub _output_ok {
	my ( $result, $outfile ) = @_;
	ok ( $result->output =~ /DONE/, "command prints 'DONE'" );
	ok ( ! $result->error, "no error in result" );
	is ( $result->run_rv, 1, "return value is 1" );
	ok ( -e $outfile, "outfile exists" );
	ok ( -s $outfile, "outfile not empty" );
}

my $infile = "$Bin/testdata/testnames.txt";
my ($fh, $outfile) = tempfile( 'CLEANUP' => 1, 'DIR' => $ENV{HOME} );

my $result = test_app( Bio::SUPERSMART::App::smrt=> [ "taxize",  "-i",   $infile, "-o", $outfile ]);
_output_ok ( $result, $outfile );

# test with wordir argument and output file name without path
my $dir = tempdir ( 'CLEANUP' => 1, 'DIR' => $ENV{HOME} );
$outfile= $dir . "/out.txt";

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "taxize",  "-i",   $infile, "-o", $outfile, "-w", $dir ]);
_output_ok ( $result, $outfile );

# test if we get error when providing incorrect infile
$infile = "notthere.txt";

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "taxize",  "-i",   $infile, "-o", $outfile ]);
ok ( $result->error, "error when infile not present" );

# test if 'Usage' is printed when providing incorrent arguments
$result = test_app( Bio::SUPERSMART::App::smrt=> [ "taxize",  "-q",   "incorrect", "-l", "argument" ]);
ok ( $result->error, "error when arguments invalid" );
ok ( $result->error =~ /Usage/, "print usage on error" );

# test root taxon expansion
$infile = "$Bin/testdata/root-taxon.txt";
$result = test_app( Bio::SUPERSMART::App::smrt=> [ "taxize",  "-i",   $infile, "-w", $dir, "-e", "species", "-o", $outfile ]);
_output_ok ( $result, $outfile );

# test with root taxon expansion given vie -r option
$result = test_app( Bio::SUPERSMART::App::smrt=> [ "taxize",  "-r",  "Homo,Felidae", "-w", $dir, "-o", $outfile ]);
_output_ok ( $result, $outfile );



