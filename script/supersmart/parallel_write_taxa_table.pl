#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Service::ParallelService 'mpi'; # can be either 'pthreads' or 'mpi';

=head1 NAME

parallel_write_taxa_table.pl - reconciles taxon labels with a taxonomy

=head1 SYNOPSYS

 $ mpi_write_taxa_table.pl --infile=<file> --root_taxon<taxon_name>  > <outfile>

=head1 DESCRIPTION

Given a file with putative taxon names, or optionally a root taxon name, attempts to link these to the pipeline's
underlying taxonomy. Writes the results as a tab-separated file to STDOUT.

=cut

# process command line arguments
my $verbosity = INFO;
my @levels = qw[varietas subspecies species genus tribe subfamily family superfamily parvorder infraorder suborder order subclass class phylum kingdom];

my ( $infile, $root_taxon );
GetOptions(
	'infile=s'   => \$infile,
	'root_taxon=s' => \$root_taxon,
	'verbose+' => \$verbosity,
	'level=s'  => \@levels,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => [qw(
		main 
		Bio::Phylo::PhyLoTA::Service::ParallelService		
	)]
    );

my $config = Bio::Phylo::PhyLoTA::Config->new;
my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;

# if 'infile' argument is given, read names from file or STDIN, clean line breaks
my @names;
if ( $infile ){
	if ( $infile eq '-' ) {
		@names = <STDIN>;
		chomp(@names);
		$log->info("read species names from STDIN");
	}
	else {
		open my $fh, '<', $infile or die $!;
		@names = <$fh>;
		chomp(@names);
		$log->info("read ".scalar(@names)." species names from $infile");
	}
}

# expand possible root taxa to taxonomic level set in ROOT_TAXA_EXPANSION 
#  or set to "Species" in the case a root taxon is given as argument, but 
#  expansion is not set in config file
my $expand_to = $config->ROOT_TAXA_EXPANSION;
	
# if a root taxon is given, add it to the taxon names and set expansion to "Species" if not set otherwise
if ( $root_taxon ) {
	push @names, $root_taxon;
	if ( ! $expand_to ){
		$expand_to = "Species"
	}
}

@names = $mts->expand_taxa( \@names , $expand_to);

# print table header to stdout in sequential mode
sequential { print join ("\t", 'name', @levels), "\n"; };

# this will take some time to do the taxonomic name resolution in the
# database and with webservices. The below code runs in parallel
my @result = pmap {
        my $name = $_; 
        
        # replace consecutive whitespaces with one
        $name =~ s/\s+/ /g ;               
        my @res = ();
        my @nodes = $mts->get_nodes_for_names($name);
        if ( @nodes ) {
                if ( @nodes > 1 ) {
                        $log->warn("found more than one taxon for name $name");
                }
                else {
                        $log->info("found exactly one taxon for name $name");
                }
                
                # for each node, fetch the IDs of all taxonomic levels of interest
                for my $node ( @nodes ) {
                        
                        # create hash of taxonomic levels so that when we walk up the
                        # taxonomy tree we can more easily check to see if we are at a
                        # level of interest
                        my %level = map { $_ => "NA" } @levels;
                        
                        # traverse up the tree
                        while ( $node ) {
                                my $tn = $node->taxon_name;
                                my $ti = $node->ti;
                                my $rank = $node->rank;
                                if ( exists $level{$rank} ) {
                                        $level{$rank} = $node->get_id;
                                } 
                                $node = $node->get_parent;
                        }                        
                        my @entry = ($name, @level{@levels});
                        push @res, \@entry;
        		}
        }
        else {
                $log->warn("couldn't resolve name $name");
        }                
        return  @res;
} @names; 


# filter taxa with duplicate species id and write the results to stdout
sequential {
	my %seen = ();
	foreach my $res ( @result ) {
		my @entry = @{$res};
		my $name =  shift @entry;
		my $ids = join "\t", @entry;
		if ( exists $seen{$ids} ) {
			$log->info("Taxon with resolved name $name already in species table. Ignoring.");
		} 
		else {
			print "$name \t $ids \n";
			$seen{$ids} = 1;	
		}
	}
}
