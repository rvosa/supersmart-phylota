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

my $config = Bio::SUPERSMART::Config->new;

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

my $matrix = "$Bin/testdata/supermatrix-pantherinae.phy";
my $taxafile = "$Bin/testdata/species-pantherinae.tsv";
my $workdir = tempdir( CLEANUP => 1, 'DIR' => $ENV{HOME} );

my ($fh, $outfile) = tempfile( 'CLEANUP' => 1, 'DIR' => $workdir );

# test inference with examl
my $result = test_app( Bio::SUPERSMART::App::smrt=> [ 'bbinfer',  '-s',   $matrix, '-t', $taxafile, '-o', $outfile, '-w', $workdir ]);
print "RESULT : " . $result->output . "\n";
_output_ok ( $result, $outfile );

# takes quite long, that's why currently commented out
# test inference with raxml from random starttree
#$result = test_app( Bio::SUPERSMART::App::smrt=> [ 'bbinfer',  '-s', $matrix, '-i', 'raxml', '-m', '-b', 10, '-o', $outfile, '-w', $workdir]);
#_output_ok ( $result, $outfile );

# test inference with exabayes
$config->EXABAYES_NUMGENS(1000);
$result = test_app( Bio::SUPERSMART::App::smrt=> [ 'bbinfer',  '-s',   $matrix, '-i', 'exabayes', '-o', $outfile, '-w', $workdir]);
_output_ok ( $result, $outfile );
