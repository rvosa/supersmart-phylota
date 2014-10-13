package Bio::Apps::Supersmart::Command::Taxize;

use strict;
use warnings;

use List::MoreUtils qw(firstidx uniq);

use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Service::ParallelService 'pthreads'; # can be either 'pthreads' or 'mpi';

use base 'Bio::Apps::GlobalCmd';
use Bio::Apps::Supersmart qw(-command);

# ABSTRACT: writes taxa table for given list of taxon names or root taxa

=head1 NAME

Taxize.pm - reconciles taxon names with information from NCBI taxonomy

=head1 SYNOPSYS

$ smrt taxize [-v ] [-w <dir>] -i <file> [-o <file>] [-e <rank>] 

=head1 DESCRIPTION

Performs taxonomic name resolution by mapping a provided list of input names onto the 
NCBI taxonomy. Optionally expands the higher input taxa to the specified lower taxon, 
e.g. to expand a named Order to its constituent species. Produces a table that lists 
all resolved (and expanded) taxa and their higher classification.

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "species.tsv";
	return (
		["infile|i=s", "file with list of taxon names", { arg => "file", mandatory => 1}],
		["outfile|o=s", "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],
		["expand_rank|e=s", "rank to which root taxa are expanded", { default => 0, arg => "rank"}],	
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	# We only have to check the 'infile' argument. 
	#  If the infile is absent or empty, abort  
	my $file = $opt->infile;
	$self->usage_error("no infile argument given") if not $file;
	$self->usage_error("file $file does not exist") unless (-e $file);
	$self->usage_error("file $file is empty") unless (-s $file);			
}

sub execute {
	my ($self, $opt, $args) = @_;
	
	my $verbosity = INFO;
	
	# collect command-line arguments
	my $infile = $opt->infile;
	my $expand_rank = $opt->expand_rank;
	my $outfile = $opt->outfile;
	(my $workdir = $opt->workdir) =~ s/\/$//g;
	$verbosity += $opt->verbose ? $opt->verbose : 0;
	
	# instantiate helper objects
	my $log = Bio::Phylo::Util::Logger->new(
		'-level' => $verbosity,
		'-class' => [ __PACKAGE__,			 
					 "Bio::Phylo::PhyLoTA::Service::ParallelService"]
    );
    
    my $config = Bio::Phylo::PhyLoTA::Config->new;
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
    
	# read names from file or STDIN, clean line breaks
	open my $fh, '<', $infile or die $!;
	my @names = <$fh>;
	chomp(@names);
	close $fh;
	$log->info("read ".scalar(@names)." species names from $infile");
	
	# expand root taxa if argument is provided
	if ( $expand_rank ) {
		@names = $mts->expand_taxa( \@names , $expand_rank);
	}
		
	# get all possible taxonomic ranks
	my @levels = reverse $mts->get_taxonomic_ranks;
		
	# this will take some time to do the taxonomic name resolution in the
	# database and with webservices. The below code runs in parallel
	my @result = pmap {
	        my $name = $_; 
	        # replace consecutive whitespaces with one
	        $name =~ s/\s+/ /g;               
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
	
	# clean the results for duplicates and rows that represent taxa with levels higher 
	# than 'Species', then write the results to the output file
	sequential {
		open my $out, '>', $workdir . '/' . $outfile or die $!;
	
		# print table header 
		print $out join ("\t", 'name', @levels), "\n";
	
		my %seen = ();
		foreach my $res ( @result ) {
			my @entry = @{$res};
			my $name =  shift @entry;
			my $ids = join "\t", @entry;
			
			# omit all taxa with higher taxonomic rank than 'Species' 
			my $highest_rank = "Species";
			my $idx = firstidx { lc ($_) eq lc ($highest_rank) } @levels;			
			my @subset = @entry[0..$idx];
			if ( uniq (@subset) == 1) {
				$log->debug("rank of taxon with resolved name $name is higher than $highest_rank, omitting");									
			}	
			# filter taxa with duplicate species id		
			elsif ( exists $seen{$ids} ) {
				$log->debug("taxon with resolved name $name already in species table, omitting");
			} 
			else {
				print $out "$name \t $ids \n";
				$seen{$ids} = 1;	
			}
		}
		close $out;		
	};
	$log->info("DONE, results written to $workdir/$outfile");
		
	return 1;    
}

1;
