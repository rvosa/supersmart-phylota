package Bio::Phylo::PhyLoTA::Service::ParallelService;
use strict;
use warnings;
use Carp 'croak';
use POSIX 'ceil';
use base 'Bio::Phylo::PhyLoTA::Service';

my $mode;
my $num_proc;

=over

=item END

MPI_Finalize need to invoked ecactly once (per mpirun execution) after all parallel 
jobs have finished, therefore this is done in the END block.

=cut

END { 
        if ( $mode eq 'mpi' ){
                MPI_Finalize();
        }
}

=item import

When modules are 'use'd, the 'import' sub is automatically called as a package method,
i.e. ParallelService->import. If the 'use' statement has additional arguments, these
are passed into the import statement. For example:

 use ParallelService 'mpi';

Results in the import sub receiving 1. package name, 2. a string determining the mode. 
The mode can either be 'mpi' or 'pthreads', if none is given, 'pthreads' is the default.

When importing this module, the subroutines 'pmap' and 'sequential' are exported to the caller's namespace, 
including thir prototypes (&@ or &). 
Pmap can be used as a drop-in replacement of the built-in perl function map. The method 'sequential' accepts a code
block which is evaluated only on the master node. Note that 'sequential' does not need to be used when in 'pthread'
mode, but it is still present for convenience. 
 
The number of processes is determined 
automatically in this class, if e.g. the script which uses the ParallelService class
is invoced with 'mpirun -np 5', the number of processes is set to 4 (4 worker nodes and a master node). In 'pthread'
mode, the number of processes is set to the number of cores present on the machine.

=cut

sub import {
	my $package = shift;
	( $mode ) = @_;	
        if ( ! $mode ) {
                $mode = 'pthreads';
                $package->logger->warn("parallel mode not provided to ParallelService. Using default mode 'pthreads'");
        }
        if ( $mode eq 'mpi' ){
                use Parallel::MPI::Simple;
                MPI_Init();       
                $num_proc = MPI_Comm_size(MPI_COMM_WORLD())-1;                
        }
        elsif ( $mode eq 'pthreads' ) {
                use threads;
                use threads::shared;
                # get and set number of processes
                use Sys::Info;
                use Sys::Info::Constants qw( :device_cpu );
                my $info = Sys::Info->new;
                my $cpu  = $info->device('CPU');
                $num_proc = $cpu->count;
        }
	my ( $caller ) = caller();
        eval "sub ${caller}::pmap(\&\@);";
	eval "*${caller}::pmap = \\&${package}::pmap;";
        
        eval "sub ${caller}::sequential(\&);";
	eval "*${caller}::sequential = \\&${package}::sequential;";
}


=item pmap_pthreads

Does a parallel version of map using pthreads. The function pmap deletages to this 
function if ParallelService was called with 'pthreads' as its argument.

=cut

sub pmap_pthreads (&@) {
	my ( $func, @data ) = @_;
	my $size = scalar @data;
        my @threads;
	my @result;
	my $inc = ceil($size/$num_proc);
	for ( my $i = 0; $i < $size; $i += $inc ) {
		my $max = ( $i + $inc - 1 ) >= $size ? $size - 1 : $i + $inc - 1;
		my @subset = @data[$i..$max];
                push @threads, threads->create( sub { map {$func->()} @subset } );
        }
	push @result, $_->join for @threads;
	return @result;
}

=item pmap_mpi

Does a parallel version of map using MPI. The function pmap deletages to this 
function if ParallelService was called with 'mpi' as its argument. Note the second
argument:

 ParallelService 'mpi' => 4;

Here, the number 4 indicates the number of WORKER nodes. Hence, the script needs 
to be called as:

 mpirun -np 5 <scriptname>

Because there is also a boss node.

=cut

sub pmap_mpi (&@) {
	my ( $func, @data ) = @_;
	my $size = scalar @data;
	my $inc = ceil($size/$num_proc);
        
	my $WORLD  = MPI_COMM_WORLD();
	my $BOSS   = 0;
	my $JOB    = 1;
	my $RESULT = 2;
	my $rank = MPI_Comm_rank($WORLD);

	# this is the boss node
	if ( $rank == 0 ) {
                # submit the data in chunks
		my $worker = 1;
		for ( my $i = 0; $i < $size; $i += $inc ) {
			my $max = ( $i + $inc - 1 ) >= $size ? $size - 1 : $i + $inc - 1;
			my @subset = @data[$i..$max];
			MPI_Send(\@subset,$worker,$JOB,$WORLD);
			$worker++;
		}
		
		# receive the result
		my @result;
                # beware of the case that there are more nodes than items to process
		for $worker ( 1 .. ceil($size/$inc) ) {
			my $subset = MPI_Recv($worker,$RESULT,$WORLD);
			push @result, @{ $subset };
		}
                return @result;
	}
	else {
	
		# receive my job and process it
		my $subset = MPI_Recv($BOSS,$JOB,$WORLD);
		my @result = map { $func->($_) } @{ $subset };
		
		# send the result back to the boss
		MPI_Send(\@result,$BOSS,$RESULT,$WORLD);
	}
}

=item sequential

This wrapper allows code to be executed in sequential mode. This is needed
if a script is called using 'mpirun' but one would like to be able to
run some code sequentially and not as many times as passed to mpirun's '-np' 
argument. Code can be simply run sequentially as follows:

my $something = sequential { some code here };

When using 'pthreads' for parallelization, the script is called
with the normal 'perl' command and thus the code in it would run sequencially
anyways unless something is sent to the ParallelService (via pmap). The function
'sequential' is therefore not needed. Thus, in mode 'pthreads', this prototype
simply evaluates the code given as argument.

=cut

sub sequential (&) {        
        my ($func) = @_;
        croak "Need mode argument!" if not $mode;
        if ( $mode eq 'pthreads' ) {
                # simply call the function
                $func->();
        }
        elsif ( $mode eq 'mpi' ) {    
                # only allow master node to perform operation
                my $WORLD  = MPI_COMM_WORLD();
                my $rank = MPI_Comm_rank($WORLD);
                if ($rank==0){
                        $func->();
                }
        }
}

=item pmap
    
This is a wrapper around whichever implementation the user has indicated. Note the
usage of 'goto'. Since the seminal work of Dutch computer scientist Edsger Dijkstra,
we all know that "go to" statements should be considered harmful
(for the original paper from 1968: http://www.u.arizona.edu/~rubinson/copyright_violations/Go_To_Considered_Harmful.html),
however, in this case it should be permissible because it is obviously labelled where
we delegate to, and it serves the purpose of erasing pmap() from the call stack, 
directly forwarding @_ to the actual handler.

=cut

sub pmap (&@) {
	croak "Need mode argument!" if not $mode;
	if ( $mode eq 'pthreads' ) {
		goto &pmap_pthreads;
	}
	elsif ( $mode eq 'mpi' ) {
		goto &pmap_mpi;
	}
}

=back

=cut

1;