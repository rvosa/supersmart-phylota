#!/usr/bin/perl
use warnings;

BEGIN {
	our @plan = $ENV{'TEST_SMRT'} ? 'no_plan' : 'skip_all' => 'env var TEST_SMRT not set';
}
use Test::More @plan;

use App::Cmd::Tester;
use FindBin '$Bin';
use File::Temp qw(tempfile);
use Bio::Phylo::IO 'parse_tree';

use Bio::SUPERSMART::App::smrt;

BEGIN { use_ok('Bio::SUPERSMART::App::smrt::Command::Classify') };

# test if run gave no errors and output file is good
sub _output_ok {
	my ( $result, $outfile ) = @_;
	ok ( $result->output =~ /DONE/, "command prints 'DONE'" );
	ok ( ! $result->error, "no error in result" );
	is ( $result->run_rv, 1, "return value is 1" );
	ok ( -e $outfile, "outfile exists" );
	ok ( -s $outfile, "outfile not empty" );
}

my @infiles = ( "$Bin/testdata/species.tsv", "$Bin/testdata/species-allranks.tsv" );

# first test if command runs and has output

foreach my $in ( @infiles ) {
	my ($fh, $outfile) = tempfile( 'CLEANUP' => 1 , 'DIR' => $ENV{HOME});
	my $result = test_app( Bio::SUPERSMART::App::smrt=> [ "classify",  "-i",   $in, "-o", $outfile ]);
	
	# check if run was ok		
	_output_ok ( $result, $outfile );
		
	# get number of taxa in the classification tree
	my $tree = parse_tree( '-file' => $outfile, '-format' => 'newick' );
	my $num_nodes = scalar ( @{$tree->get_terminals} ); 
	cmp_ok ( $num_nodes, ">", 0, "number of nodes in tree not zero" );			
}

# test if we get error when providing incorrect infile
my $infile = "notthere.txt";

my $result = test_app( Bio::SUPERSMART::App::smrt=> [ "classify",  "-i",   $infile ]);
ok ( $result->error, "error when infile not present" );

# test if 'Usage' is printed when providing incorrent arguments
$result = test_app( Bio::SUPERSMART::App::smrt=> [ "classify",  "-q",   "incorrect", "-k", "argument" ]);
ok ( $result->error, "error when arguments invalid" );
ok ( $result->error =~ /Usage/, "print usage on error" );

