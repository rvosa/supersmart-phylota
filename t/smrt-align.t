#!/usr/bin/perl
use warnings;

use Test::More 'no_plan';
use App::Cmd::Tester;
use FindBin '$Bin';
use File::Temp qw(tempfile tempdir);
use Bio::SUPERSMART::App::smrt;

BEGIN {
	use Bio::Phylo::Util::Logger ':levels';
	our $log = Bio::Phylo::Util::Logger->new( 
		'-level' => INFO, 
		'-class' => [
			'main',
			'Bio::Phylo::PhyLoTA::Service::ParallelService',
		]
	);
}

BEGIN { use_ok('Bio::SUPERSMART::App::smrt::Command::Align') };

my $workdir = tempdir();
my $outfile = $workdir . "/aligned.txt";
my $infile  = "$Bin/testdata/species-fishes.tsv";

$log->info("Workdir: $workdir");

my @result = `$Bin/../script/smrt align -w $workdir -i $infile -o $outfile`;

ok( $? == 0, "no error code" );
ok( $result[-1] =~ /DONE/, "command prints DONE" );
ok( -e $outfile, "$outfile exists" );
ok( -s $outfile, "$outfile not empty" );
