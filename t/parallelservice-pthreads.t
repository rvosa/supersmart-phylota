use strict;
use warnings;

use Test::More 'no_plan';


BEGIN { use_ok('Bio::Phylo::PhyLoTA::Service::ParallelService', 'pthreads') };

cmp_ok ( num_workers(), ">", 0, "More than one worker in default configuration");

my @letters;
sequential{ @letters = ('A'..'K') };
my @result = pmap {sleep(1); print "Processing item $_ \n"; return $_++; } @letters ;
sequential{ cmp_ok (scalar @result, '==', 11, "Returned array has correct length") };
	

# Do two calculations: One wih four nodes and one with one node. 
#  The parallel calculation should be faster!

num_workers(1);
my $before = time();
@result = pmap {sleep(1); print "Processing item $_ \n"; return $_++; } @letters ;
my $after = time();
my $elapsed_sequential = $after - $before;

num_workers(4);
$before = time();
@result = pmap {sleep(1); print "Processing item $_ \n"; return $_++; } @letters ;
$after = time();

my $elapsed_parallel = $after - $before;
	
cmp_ok ( $elapsed_parallel, '<', $elapsed_sequential, "parallel calculation faster than sequential");
	

