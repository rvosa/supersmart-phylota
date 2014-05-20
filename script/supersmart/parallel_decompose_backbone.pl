#!/usr/bin/perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Data::Dumper;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Service::ParallelService 'mpi'; # can be either 'pthreads' or 'mpi';

=head1 NAME

decompose_backbone.pl - decomposes backbone in monophyletic groups, writes alignments

=head1 SYNOPSYS

 $ perl decompose_backbone.pl -t <taxa> -l <list> -b <tree> -w <dir> [--verbose]

=head1 DESCRIPTION

Traverses the backbone C<tree> to find the nearest monophyletic clade that groups
the exemplar leaves. In the default case, the clade is the genus that subtends the
two exemplars provided they are monophyletic in the backbone tree. If they are not, we
traverse upwards to find the nearest monophyletic set of genera.

Each clade is then expanded into its constituent set of species on the basis of the
taxon mapping file C<taxa>. For those sets of species, the C<list> of alignment values
is evaluated, and for each alignment whose average divergence does not exceed 
CLADE_MAX_DISTANCE but whose density of species in the sets does exceed CLADE_MIN_DENSITY
the sequences for the focal species are written to a new file, in a directory that
groups them with the other relevant alignments for that clade.

=cut

# process command line arguments
my $verbosity = WARN;
my ( $taxa, $list, $backbone, $workdir );
GetOptions(
	'taxa=s'     => \$taxa,
	'list=s'     => \$list,
	'backbone=s' => \$backbone,
	'workdir=s'  => \$workdir,
	'verbose+'   => \$verbosity,
);

# instantiate helper objects
my $mt = 'Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa';
my $ts = 'Bio::Phylo::PhyLoTA::Service::TreeService';
my $config = Bio::Phylo::PhyLoTA::Config->new;
my $logger = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => ['main' ,
                     'Bio::Phylo::PhyLoTA::Service::ParallelService',
        ], 
);

# parse backbone tree
$logger->info("going to read backbone tree $backbone");
my $tree = parse_tree(
	'-format'     => 'newick',
	'-file'       => $backbone,
	'-as_project' => 1,
);

# parse taxon mapping
$logger->info("going to read taxa mapping $taxa");
my @records = $mt->parse_taxa_file($taxa);


# decompose tree into clades and get the sets of species
my @set = $ts->_make_clade_species_sets($tree, @records);

# now read the list of alignments
my @alignments;
sequential {
	$logger->info("going to read list of alignments $list");
	open my $fh, '<', $list or die $!;
	while(<$fh>) {
		chomp;
		push @alignments, $_ if /\S/ and -e $_;
	}
};

# write suitable alignments to their respective clade folders


pmap {
        my $aln = $_;
	$logger->debug("assessing $aln");
	my %fasta = $mt->parse_fasta_file($aln);
	my $dist  = $mt->calc_mean_distance(%fasta);
	my $nchar = $mt->get_nchar(%fasta);
	my $mdens = $config->CLADE_MIN_DENSITY;
	if ( $dist <= $config->CLADE_MAX_DISTANCE ) {
	
		# assess for each set whether we have enough density
		for my $i ( 0 .. $#set ) {
			my @species = @{ $set[$i] };
			my %seq;
			my $distinct = 0;
			for my $s ( @species ) {
				my @s = grep { /taxon\|$s[^\d]/ } keys %fasta;
				if ( @s ) {
					$seq{$s} = \@s;
					$distinct++;
				}
			}
			
			# the fraction of distinct, sequenced species is high enough,
			# and the total number exceeds two (i.e. there is some topology to resolve)
			if ( ($distinct/scalar(@species)) >= $mdens && $distinct > 2 ) {
				$logger->info("$aln is informative and dense enough for clade $i");
				my ( $fh, $seed ) = make_handle( $i, $aln );
				for my $s ( @species ) {
					if ( $seq{$s} ) {
						for my $j ( 0 .. $#{ $seq{$s} } ) {
							my $seq = $seq{$s}->[$j];
							print $fh '>', $seq, "\n", $fasta{$seq}, "\n";
						}
					}
				}
			}
		}
	}
	else {
		$logger->info("$aln is too divergent: $dist > ".$config->CLADE_MAX_DISTANCE);
	}
} @alignments;

sub make_handle {
	my ( $i, $aln ) = @_;
	if ( not -d "$workdir/clade$i" ) {
		mkdir "$workdir/clade$i";
	}
	my ( $volume, $dir, $base ) = File::Spec->splitpath($aln);
	open my $fh, '>', "$workdir/clade$i/$base" or die $!;
	$base =~ s/\.fa$//;
	return $fh, $base;
}
