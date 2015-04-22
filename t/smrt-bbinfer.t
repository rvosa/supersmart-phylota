#!/usr/bin/perl
use warnings;

use Test::More 'no_plan';
use App::Cmd::Tester;
use FindBin '$Bin';
use File::Temp qw(tempfile tempdir);

use Bio::SUPERSMART::App::smrt;

my $config = Bio::Phylo::PhyLoTA::Config->new;

BEGIN { use_ok('Bio::SUPERSMART::App::smrt::Command::BBinfer') };

# test if run gave no errors and output file is good
sub _output_ok {
	my ( $result, $outfile ) = @_;
	print "ERROR : ";
	print $result->error . "\n";
	#die;
	#ok ( $result->output =~ /DONE/, "command prints 'DONE'" );
	ok ( ! $result->error, "no error in result" );
	is ( $result->run_rv, 1, "return value is 1" );
	ok ( -e $outfile, "outfile exists" );
	ok ( -s $outfile, "outfile not empty" );
}

my $matrix = "$Bin/testdata/supermatrix.phy";

my $workdir = tempdir( CLEANUP => 1, 'DIR' => $ENV{HOME} );

my ($fh, $outfile) = tempfile( 'CLEANUP' => 1, 'DIR' => $workdir );

# test inference with examl and bootstrapping
my $result = test_app( Bio::SUPERSMART::App::smrt=> [ 'bbinfer',  '-s',   $matrix, '-o', $outfile, '-w', $workdir ]);
_output_ok ( $result, $outfile );

# test inference with raxml 
$result = test_app( Bio::SUPERSMART::App::smrt=> [ 'bbinfer',  '-s', $matrix, '-i', 'raxml', '-o', $outfile, '-w', $workdir]);
_output_ok ( $result, $outfile );

# test inference with exabayes
$config->EXABAYES_NUMGENS(1000);
$result = test_app( Bio::SUPERSMART::App::smrt=> [ 'bbinfer',  '-s',   $matrix, '-i', 'exabayes', '-o', $outfile, '-w', $workdir]);
_output_ok ( $result, $outfile );
