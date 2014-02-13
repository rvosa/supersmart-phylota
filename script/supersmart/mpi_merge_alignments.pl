#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use File::Temp 'tempfile';
use Getopt::Long;
use Parallel::MPI::Simple;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use constant HEAD_NODE => 0;
use constant ALIGNMENT_FILES_SUBSET => 1;
use constant PROTID_FOR_SEED_GI => 2;
use constant CLUSTER_SUBSET => 3;
use constant PROFILE_ALIGNMENTS => 4;

=head1 NAME

mpi_merge_alignments.pl - merges alignments by various means

=head1 SYNOPSYS

 $ mpi_merge_alignments.pl --list=<file> --stem=<base_name> > <outfile>

=head1 DESCRIPTION

Given a file containing a list of alignments, merges these alignments into orthologous
sets. Writes the names of the newly created, merged alignment files to STDOUT.

=cut

# process command line arguments
my ( $list, $stem );
my $verbosity = WARN;
GetOptions(
	'list=s'   => \$list,
	'stem=s'   => \$stem,
	'verbose+' => \$verbosity,
);

# instantiate helper objects
Bio::Phylo::Util::Logger->VERBOSE(
	'-class' => [
		'main',
		'Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector',
		'Bio::Phylo::PhyLoTA::Service::SequenceGetter',
	],
	'-level' => $verbosity,
);
my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
my $sg  = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;

# here we start in parallel mode
MPI_Init();
my $rank = MPI_Comm_rank(MPI_COMM_WORLD);

# this is executed by the boss node
if ( $rank == 0 ) {

    # read list of files
    my @list = read_file($list);

    # MPI_Comm_size returns the total number of nodes. Because we have one
    # head node this needs to be - 1
    my $nworkers = MPI_Comm_size(MPI_COMM_WORLD) - 1;
    INFO "we have $nworkers nodes available";
    
    # first look up the InParanoid protein ID for the seed GIs
    my @result = dispatch_job( $nworkers, ALIGNMENT_FILES_SUBSET, PROTID_FOR_SEED_GI, @list );
    
    # now we have a list where each element is a hash ref whose keys are seed GIs, values are
    # InParanoid protein IDs. We need to do an index inversion so that we get
    # all seed GIs for a given protein ID.
    my %seed_gis_for_protid;
    my %orthologs_for_protid;
    
    # @result is an array of hash references. there are as many hash references as there are
    # available workers
    for my $hash ( @result ) {
        for my $seed_gi ( keys %{ $hash } ) {
            my $protid = $hash->{$seed_gi};
            
            # here we do the index inversion
            $seed_gis_for_protid{$protid} = [] if not $seed_gis_for_protid{$protid};
            push @{ $seed_gis_for_protid{$protid} }, $seed_gi;
            
            # for all protein IDs we also store a lookup of all their orthologs, for later use
            if ( not $orthologs_for_protid{$protid} ) {
                INFO "looking up orthologs for protein ID $protid";
                my %orthologs = map { $_ => 1 } $sg->get_orthologs_for_protein_id($protid);
                $orthologs_for_protid{$protid} = \%orthologs;
            }
        }
    }
    
    # here we cluster the protein IDs by orthology, @clusters is a list of lists of
    # clustered protein IDs
    my @protids = keys %seed_gis_for_protid;
    my @clusters = cluster_protein_ids(\%orthologs_for_protid, @protids);
    
    # translate each cluster back to file names
    my @translated;
    for my $c ( @clusters ) {
        my %cluster;
        for my $pid ( @{ $c } ) {
            my @files = @{ $seed_gis_for_protid{$pid} };
            $cluster{$_} = 1 for @files;
        }
        push @translated, [ keys %cluster ];
    }
    
    # profile align the files
    my @filelists = dispatch_job( $nworkers, CLUSTER_SUBSET, PROFILE_ALIGNMENTS, @translated );
    for my $list ( @filelists ) {
        for my $file ( @{ $list } ) {
            print $file, "\n";
        }
    }
}

# this is executed by worker nodes
else {
    
    # lookup InParanoid protein IDs for seed GIs
    my $subset = MPI_Recv(0,ALIGNMENT_FILES_SUBSET,MPI_COMM_WORLD);
    INFO "worker $rank has received ".scalar(@$subset)." files to analyse";
    my %result;
    for my $file ( @{ $subset } ) {
        if ( $file =~ /(\d+)\.fa/ ) {
            my $seed_gi = $1;
            if ( my $protid = $sg->get_protid_for_seed_gi($seed_gi) ) {
                $result{$file} = $protid;
                INFO "worker $rank found protein ID $protid for seed GI $seed_gi";
            }
            else {
                WARN "worker $rank found no protein ID for seed GI $seed_gi";
            }
        }
        else {
            WARN "worker $rank couldn't parse seed GI from file name $file";
        }
    }
    MPI_Send(\%result,HEAD_NODE,PROTID_FOR_SEED_GI,MPI_COMM_WORLD);

    # perform recursive profile alignment
    my $clusters = MPI_Recv(0,CLUSTER_SUBSET,MPI_COMM_WORLD);
    INFO "worker $rank has received ".scalar(@$clusters)." clusters to align";
    my @result;
    my $counter = 1;
    for my $cluster ( @{ $clusters } ) {
        INFO "worker $rank will align cluster $counter";
        my @files = @{ $cluster };
        my $outfile = "$stem-$rank-$counter.fa";
        copy(shift @files, $outfile);
        if ( @files ) {
	    for my $i ( 0 .. $#files ) {
                INFO "worker $rank is aligning file ".($i+2)." of cluster $counter";
                my $result = $sg->profile_align_files($outfile,$files[$i]);
                open my $fh, '>', $outfile or die $!;
                print $fh $result;
                close $fh;
            }
        }
        push @result, $outfile;
        $counter++;
    }
    MPI_Send(\@result,HEAD_NODE,PROFILE_ALIGNMENTS,MPI_COMM_WORLD);
}

# we need to exit cleanly by making this call on all nodes
MPI_Finalize();

=begin comment

Clusters a list of protein IDs by orthology. Returns a list of lists
of clustered protein IDs

=cut comment

sub cluster_protein_ids {
    my ( $orthologs_for_protid, @protids ) = @_;
    my ( %seen, @clusters );
    for my $i ( 0 .. $#protids ) {
        my $protid1 = $protids[$i];
        
        # protein IDs are no longer a candidate for clustering
        # if they've already been clustered, so one $i gets higher
        # in the list there are things to skip here, because they've
        # already been clustered in an earlier inner loop
        if ( not $seen{$protid1} ) {
            
            # initialize focal cluster
            my @cluster = ( $protid1 );
            
            # inner loop, here we compare all remaining protein IDs
            # with the focal ID
            for my $j ( ( $i + 1 ) .. $#protids ) {
                my $protid2 = $protids[$j];
                
                # some of the remaining protein IDs may have already
                # been clustered, we skip over these
                if ( not $seen{$protid2} ) {
                    
                    # the values in this hash are *all* protein IDs
                    # that are orthologous with the focal one. We
                    # check to see if $protid2 is among them
                    if ( $orthologs_for_protid->{$protid1}->{$protid2} ) {
                        push @cluster, $protid2;
                        $seen{$protid2}++;
                    }
                }
            }
            push @clusters, \@cluster;
            $seen{$protid1}++;
        }
    }
    return @clusters;
}

=begin comment

Dispatches an MPI job

=cut comment

sub dispatch_job {
    my ( $nworkers, $send_tag, $return_tag, @list ) = @_;
    
    # for each worker we dispatch an approximately equal chunk of list elements.
    my $nfiles = int( scalar(@list) / $nworkers );
    INFO "each worker will process about $nfiles elements";
    
    # this is the starting index of the subset
    my $start = 0;
    
    # iterate over workers to dispatch jobs
    for my $worker ( 1 .. $nworkers ) {
	
	# the ending index of the subset is either the starting index + chunk size or the
	# highest index in the array, whichever fits
	my $end  = ( $start + $nfiles ) > $#list ? $#list : ( $start + $nfiles );
	
	# create subset
	my @subset = @list[ $start .. $end ];
	INFO "dispatching ".scalar(@subset)." names (index: $start .. $end) to worker $worker";
	MPI_Send(\@subset,$worker,$send_tag,MPI_COMM_WORLD);
	
	# increment starting index
	$start += $nfiles + 1;
    }
    
    # iterate over workers to receive results
    my @results;
    for my $worker ( 1 .. $nworkers ) {
	
	# get the result
	my $result = MPI_Recv($worker,$return_tag,MPI_COMM_WORLD);
        push @results, $result;
	INFO "received results from worker $worker";	        
    }
    return @results;
}

=begin comment

Slurps a file into an array

=cut comment

sub read_file {
    my $list = shift;
    my $fh; # file handle
    
    # may also read from STDIN, this so that we can pipe 
    if ( $list eq '-' ) {
            $fh = \*STDIN;
            DEBUG "going to read file names from STDIN";
    }
    else {
            open $fh, '<', $list or die $!;
            DEBUG "going to read file names from $list";
    }
    
    # read lines into array
    my @list = <$fh>;
    chomp @list;
    
    # return result
    return @list;
}