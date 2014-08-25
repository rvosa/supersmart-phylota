use strict;
use warnings;

use FindBin '$Bin';
use Test::More 'no_plan';

(my $scriptname =  $0) =~ s/^.+\///g;

# do nothing if script was not called with 'mpirun'
if ( $ARGV[0] ){
	if ( $ARGV[0] eq "-p"){
		
		BEGIN { use_ok('Bio::Phylo::PhyLoTA::Service::ParallelService', 'mpi') };
		my @numbers;
		sequential{ @numbers = (1..10) };
		my @result = pmap { print "Processing item $_\n"; return $_++; } @numbers ;
		sequential{ cmp_ok (scalar @result, '==', 10, "Returned array has correct length") };
	}
}			
else {
	# call script with mpirun and flag for parallel mode
	system ( "mpirun -np 4 perl $Bin/$scriptname -p");	
}	
