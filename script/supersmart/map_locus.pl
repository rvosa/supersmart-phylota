#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;

# process command line arguments
my $verbosity = WARN;
my $tree_file;      # newick tree file with bootstrap values
my $alignments;     # list of alignments
my $locus = 'ITS1'; # locus name
my $bootmin = 0.85; # minimum bootstrap value
my $sizemin = 10;   # minimum clade size
my $cover   = 0.75; # minimum cover of break node
my $workdir;        # directory to write break node alns to
my $taxatable;      # species-to-genus mapping
my $help; # this script is too complicated to remember
GetOptions(
	'verbose+'     => \$verbosity,
	'tree=s'       => \$tree_file,
	'alignments=s' => \$alignments,
	'locus=s'      => \$locus,
	'bootmin=f'    => \$bootmin,
	'sizemin=i'    => \$sizemin,
	'cover=f'      => \$cover,
	'workdir=s'    => \$workdir,
	'taxatable=s'  => \$taxatable,
	'help|?'       => \$help,
);

# provide help message
if ( $help ) {
print <<"HELP";
Usage: $0 -tree <bootstrap newick tree> -a <list of alignments>
          -workdir <dir name> -taxa <ncbi taxa mapping> [-l <locus name>] 
          [-b <min bootstrap>] [-s <min clade size>] [-c <min taxon cover>]
          > <labeled tree>
HELP
exit 0;
}

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# sanity checks
$log->info("ARGS: going to concatenate $locus alignments");
$log->info("ARGS: going to read tree from $tree_file");
$log->info("ARGS: going to read alignments list from $alignments");
$log->info("ARGS: minimum bootstrap value to break on is $bootmin");
$log->info("ARGS: minimum clade size to break on is $sizemin");
$log->info("ARGS: minimum acceptable break node coverage is $cover");
$log->info("ARGS: will write concatenated alignment files in $workdir");

# instantiate helper objects
my $tree = parse_tree(
	'-format'     => 'newick',
	'-file'       => $tree_file,
	'-as_project' => 1,
);
my $sg  = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
my $mts = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;

# one-to-many taxon to file mapping, created going through the list
# of alignments and extracting which taxa occur in each FASTA file
my %files_for;
{
	$log->info("going to read alignment list from $alignments");
	open my $fh, '<', $alignments or die $!;
	while(my $aln = <$fh>) {
		chomp $aln;
		if ( $aln =~ /\.$locus\./ ) {
			$log->debug("file $aln contains $locus data");
			open my $alnfh, '<', $aln or die $!;
			while(<$alnfh>) {
				if ( /taxon\|(\d+)/ ) {
					my $taxon = $1;
					$files_for{$taxon} = [] unless $files_for{$taxon};
					push @{ $files_for{$taxon} }, $aln;
					$log->debug("taxon $taxon occurs in file $aln");
				}
			} 
		} 
	}
	$log->info("done reading alignments list");
}

# one-to-one species to genus mapping, created by reading the standard 
# TNRS mapping output file
my %genus_for;
{
    $log->info("going to read taxa mapping from $taxatable");
    open my $fh, '<', $taxatable or die $!;
    my @header;
    LINE: while(<$fh>) {
        chomp;
        my @line = split /\t/, $_;
        if ( not @header ) {
            @header = @line;
            next LINE;
        }
        my ( $species, $genus );
        for my $i ( 0 .. $#header ) {
            if ( $header[$i] eq 'species' ) {
                $species = $line[$i];
            }
            if ( $header[$i] eq 'genus' ) {
                $genus = $line[$i];
            }
        }
        $genus_for{$species} = $genus;
    }
    $log->info("done reading taxa mapping");
}

# traverse the tree, write out concatenated alignments for each node
$tree->visit_depth_first(
	'-pre' => sub {
		my $node = shift;

		# if node is a tip, store all files that contain that
		# tip and start assembling a growing list of descendants
		# for the interior nodes
		if ( $node->is_terminal ) {
			my $id = $node->get_name;
			my $files = $files_for{$id} || [];
			my %files = map { $_ => [ $id ] } @{ $files };
			$node->set_generic( 'files' => \%files );
			$node->set_generic( 'tips'  => [ $id ] );

			# not sure how to proceed in the face of incomplete data.
			# we now accept incomplete coverage, though we should
			# *at least* be retrieving all tips from the breakable
			# tree.
			if ( not $files_for{$id} or not scalar @{ $files_for{$id} } ) {
				$log->warn("no $locus files for ID $id");
			}
		}
	},
	'-post' => sub {
		my $node = shift;
		my $id = $node->get_internal_name;
		my @children = @{ $node->get_children };

		# node is internal
		if ( @children ) {
			my %merged_files;
			my @merged_tips;

			# iterate over immediate children
			for my $child ( @children ) {

				# grow the list of descendants
				push @merged_tips, @{ $child->get_generic('tips') };

				# get the mapping of files that contain descendants
				# of current child, merge it with mappings for
				# for other children 
				my $child_files = $child->get_generic('files');
				for my $child_file ( keys %{ $child_files } ) {
					$merged_files{$child_file} = [] if not $merged_files{$child_file};
					my @ids = @{ $child_files->{$child_file} };
					push @{ $merged_files{$child_file} }, @ids;
				}
			}

			# decide if this is a node that becomes a subalignment
			my $boot = $node->get_branch_length || 1; # set to 1 if root
			if ( $boot >= $bootmin && scalar @merged_tips >= $sizemin ) {
				$log->info("node $id has more than $sizemin descendants and bootstrap is $boot");

				# we break the tree at this point if at least
				# all the descendants occur in one of the PI files, in
				# which case we remove all the used files

				# select the PI files, i.e. more than 3 taxa
				my @candidates = grep { scalar @{ $merged_files{$_} } > 3 } keys %merged_files;

				# get the tips that these files cover
				my %covered = map { $_ => 1 } map { @{ $merged_files{$_} } } @candidates;
				my $fraction = scalar(keys %covered) / scalar(@merged_tips);
				if ( $fraction > $cover ) {
					$log->info("enough coverage (${fraction} > ${cover}), will write alignments");

					# write fasta matrix
					my $outfile = write_concatenated(
						'locus' => $locus,
						'id'    => $id, # internal_name
						'files' => \@candidates,
						'tips'  => \@merged_tips,
					);

					# don't write the same alignment twice!
                                        delete @merged_files{@candidates};

					# remove seen tips
                                        @merged_tips = grep { ! $covered{$_} } @merged_tips;
					
					# this so that we can later splice the subtrees into the
					# big one, by gluing the newick string of the small tree
					# into the large one
					$node->set_name( $id );
                                        $log->info("wrote alignments for $id to $outfile");  
				}
				else {
					$log->info("only $fraction tips are present in files");
				}
			}

                        # store the grown sets
                        $node->set_generic( 'files' => \%merged_files );
                        $node->set_generic( 'tips'  => \@merged_tips  );
		}
	},
);

# capture this so we can paste subtrees onto labeled nodes
print $tree->to_newick( '-nodelabels' => 1 );

sub write_concatenated {
	my %args = @_;

	# this way we only retain the species in the alignments that belong to genera
	# subtended by the focal breakable node
	my %genera_to_keep = map { $_ => 1 } map { $genus_for{$_} }  @{ $args{'tips'} };

	# iterate over files
	my %matrix;
	for my $file ( @{ $args{'files'} } ) {
		$log->info("including file $file in $locus-$args{id}");

		# parse fasta file, at this point it might contain tips
		# that are outside of the current node in the breakable tree
		# (we call this "polyphyletic" further up in this file), as
		# well as tips that should be subtended by the current node
		# but aren't because they weren't selected as exemplars in
		# the creation of the breakable tree. The first category of
		# tips we want to remove, the second we want to keep.
		my %fasta = $mts->parse_fasta_file( $file );

		# get coordinates where $locus's seed GI occurs
                my $gi;
                if ( $file =~ /\/(\d+)[^\/]+$/ ) {
                        $gi = $1;
			$log->info("seed GI is $gi");
                }
                my ( $seq ) = map { $fasta{$_} } grep { /gi\|$gi\|/ } keys %fasta;
                my @indices = $sg->get_aligned_locus_indices( $gi, $locus, $seq, 'verify' );

                # extract subsequences
		my @delete;
                DEFLINE: for my $defline ( keys %fasta ) {

                        # delete taxa that occur paraphyletically in the breakable tree
                        if ( $defline =~ /taxon\|(\d+)/ ) {
                                my $taxon = $1;

                                # only keep taxa that belong to the genera to which the exemplar
                                # species subtended by the focal node belong as well
				if ( not $genera_to_keep{ $genus_for{$taxon} } ) {
                                	push @delete, $defline;
					next DEFLINE;
				}
                        }

                        my @seq = split //, $fasta{$defline};
                        my @subseq = @seq[@indices];
                        $fasta{$defline} = join '', @subseq;
                }
		delete @fasta{@delete};
		%fasta = degap(%fasta);

		# calculate number of characters in degapped matrix
		my ($nchar) = map { length $_ } values %fasta;
		$log->info("alignment $file now has $nchar characters");

		# calculate number of characters in combined matrix
		my ($totnchar) = map { length $_ } values %matrix;

		# for any taxa that we now see for the first time in this
		# degapped matrix we need to grow the total alignment until
		# they are flush with the rest
		if ( defined $totnchar and $totnchar ) {
			$matrix{$_} .= '?' x $totnchar for grep { ! $matrix{$_} } keys %fasta;
                	$log->info("total matrix now has $totnchar characters");
		}

		# concatenate the degapped matrix with what we have so far
		$matrix{$_} .= $fasta{$_} for keys %fasta;

		# add missing to any rows that weren't in the degapped matrix
		$matrix{$_} .= '?' x $nchar for grep { ! $fasta{$_} } keys %matrix;
	}

	# write outfile
	my $outfile = "${workdir}/${locus}-$args{id}.fa";
	$log->info("going to write concatenated alignments to $outfile");
	open my $fh, '>', $outfile or die $!;
	for my $defline ( keys %matrix ) {
		print $fh '>', $defline, "\n", $matrix{$defline}, "\n";
	}
	close $fh;
	return $outfile;
}

sub degap {
	$log->info("removing all-gap columns");
        my %fasta  = @_;
        my @matrix = map { [ split //, $_ ] } values %fasta;
        my $nchar = $#{ $matrix[0] };
        my %indices;
        SITE: for my $i ( 0 .. $nchar ) {
                for my $row ( @matrix ) {
                        $indices{$i} = 1;

			# we skip to the next column if this cell is a true
			# state. consequently we won't end up deleting this
			# index, which means it is retained in the final
			# matrix
                        next SITE if $row->[$i] ne '-' and $row->[$i] ne '?';
                }
                delete $indices{$i};
        } 
        my @keep = sort { $a <=> $b } keys %indices;
        for my $defline ( keys %fasta ) {
                my @seq = split //, $fasta{$defline};
                my @subseq = @seq[@keep];
                $fasta{$defline} = join '', @subseq;
        }
        return %fasta;
}
