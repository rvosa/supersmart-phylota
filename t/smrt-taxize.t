#!/usr/bin/perl
use warnings;

use Test::More 'no_plan';
use App::Cmd::Tester;
use FindBin '$Bin';
use File::Temp qw(tempfile tempdir);

use Bio::SUPERSMART::App::smrt;

my $infile = "$Bin/testdata/testnames.txt";
my ($fh, $outfile) = tempfile( 'CLEANUP' => 1 );

my $result = test_app( Bio::SUPERSMART::App::smrt=> [ "taxize",  "-i",   $infile, "-o", $outfile ]);

ok ( $result->output =~ /DONE/, "command prints 'DONE'" );
ok ( ! $result->error, "no error in result" );
is ( $result->run_rv, 1, "return value is 1" );
ok ( -e $outfile, "outfile exists" );
ok ( -e $outfile, "outfile not empty" );

# test with wordir argument and output file name without path
my $dir = tempdir ( 'CLEANUP' => 1 );
$outfile= "out.txt";

$result = test_app( Bio::SUPERSMART::App::smrt=> [ "taxize",  "-i",   $infile, "-o", $outfile, "-w", $dir ]);

ok ( $result->output =~ /DONE/, "command prints 'DONE'" );
ok ( ! $result->error, "no error in result" );
is ( $result->run_rv, 1, "return value is 1" );
ok ( -e $outfile, "outfile exists" );
ok ( -e $outfile, "outfile not empty" );

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
$result = test_app( Bio::SUPERSMART::App::smrt=> [ "taxize",  "-i",   $infile, "-w", $dir, "-e", "species" ]);
ok ( $result->output =~ /DONE/, "command prints 'DONE'" );
ok ( ! $result->error, "no error in result" );
is ( $result->run_rv, 1, "return value is 1" );
ok ( -e $outfile, "outfile exists" );
ok ( -e $outfile, "outfile not empty" );
