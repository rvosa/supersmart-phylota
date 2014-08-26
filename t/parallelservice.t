use strict;
use warnings;

use FindBin '$Bin';
use Test::More 'no_plan';

(my $scriptname =  $0) =~ s/^.+\///g;

# do nothing if script was not called with 'mpirun'
if ( $ARGV[0] ){
	if ( $ARGV[0] eq "-p"){
		
		BEGIN { use_ok('Bio::Phylo::PhyLoTA::Service::ParallelService', 'mpi') };
		my @letters;
		sequential{ @letters = ('A'..'Z') };
		my @result = pmap { print "Processing item $_ \n"; return $_++; } @letters ;
		sequential{ cmp_ok (scalar @result, '==', 26, "Returned array has correct length") };
	}
}			
else {
	# call script with mpirun and flag for parallel mode
	system ( "mpirun -np 4 perl $Bin/$scriptname -p");	
}	
