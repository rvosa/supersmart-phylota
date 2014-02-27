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
my $config = Bio::Phylo::PhyLoTA::Config->new;
my $logger = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# parse backbone tree. XXX note that this tree must
# be rooted, while the output of examl is not.
$logger->info("going to read backbone tree $backbone");
my $tree = parse_tree(
	'-format'     => 'newick',
	'-file'       => $backbone,
	'-as_project' => 1,
);

# parse taxon mapping
$logger->info("going to read taxa mapping $taxa");
my @records = $mt->parse_taxa_file($taxa);

# traverse tree, retrieve monophyletic genera, which we will compare to all
$logger->info("going to identify putatively monophyletic genera");
my ( %monophyletic, %all );
$tree->visit_depth_first(
	'-post' => sub {
		my $node = shift;
		
		# for the tips, we keep a running tally of all seen genera. if a genus 
		# occurs on more than one tip we need to know if the genus is monophyletic
		if ( $node->is_terminal ) {
			my $id = $node->get_name;
			my ($genus) = map { $_->{'genus'} } grep { $_->{'species'} == $id } @records;
			$node->set_generic( 'species' => [ $id ], 'genera' => { $genus => 1 } );
			
			# in the end, this will be either 1 (for monotypic genera), or 2 (for
			# exemplar genera).
			$all{$genus}++;
		}
		else {
		
			# here we simply carry over the species and genera from the
			# children to the focal node. in the second pass we will use
			# this to find the shallowest node that subtends all species
			# of a paraphyletic genus
			my ( @s, %g );						
			for my $child ( @{ $node->get_children } ) {
				push @s, @{ $child->get_generic('species') };
				my %genus = %{ $child->get_generic('genera') };
				$g{$_} += $genus{$_} for keys %genus;
			}
			$node->set_generic( 'species' => \@s, 'genera' => \%g );
			
			# the node subtends two species that all belong to the same genus. 
			# hence, the node is monophyletic.
			if ( scalar(@s) >= 2 and scalar keys %g == 1 ) {
				$monophyletic{$_} += $g{$_} for keys %g;
			}
		}
	}
);

# all the genera that aren't putatively monophyletic but that do have two
# members are therefore paraphyletic. in the second pass we lump all the
# putatively monophyletic genera that nest inside paraphyletic ones within
# the mixed set of mono/para.
$logger->info("going to reconstruct paraphyly");
my %paraphyletic = map { $_ => 1 } grep { !$monophyletic{$_} && $all{$_} == 2 } keys %all;
my @genera;
$tree->visit_depth_first(
	'-post' => sub {
		my $node = shift;
		if ( $node->is_internal ) {
			my %g = %{ $node->get_generic('genera') };
			
			# the node is a paraphyletic mrca if it is the shallowest
			# node where a paraphyletic genus occurs and where it
			# subtends the two exemplars from that genus
			my $is_para;
			for my $genus ( keys %g ) {
				if ( $paraphyletic{$genus} and $g{$genus} == 2 ) {
					$is_para++;
				}
			}
			
			# if the node is paraphyletic we store ALL its subtended
			# genera, removing them from the set of putative monophyletic
			# genera as well as from the paraphyletic ones. we remove
			# from the former set so that after this traversal we don't
			# put monophyletic genera that nest inside paraphyletic ones
			# in a separate set, and we remove the latter set so that
			# deeper nodes that also subtend the paraphyletic genera
			# don't trigger processing.
			if ( $is_para ) {
				my @g = keys %g;
				push @genera, \@g;
				delete @paraphyletic{@g};
				delete @monophyletic{@g};
			}	
		}	
	}
);
# all remaining ones become their own set
push @genera, [ $_ ] for keys %monophyletic;

# now resolve the nesting of paraphyletic genera
my %index;
for my $i ( 0 .. $#genera ) {
	$index{$_} = $i for @{ $genera[$i] };
}

# make species sets
my @set;
for my $i ( keys %{ { map { $_ => 1 } values %index } } ) {
	my @g = @{ $genera[$i] };
	my @s = map { $mt->get_species_for_taxon( 'genus' => $_, @records ) } @g;
	push @set, \@s if scalar(@s) > 2;
}

# now read the list of alignments
my @alignments;
{
	$logger->info("going to read list of alignments $list");
	open my $fh, '<', $list or die $!;
	while(<$fh>) {
		chomp;
		push @alignments, $_ if /\S/ and -e $_;
	}
}

# write suitable alignments to their respective clade folders
for my $aln ( @alignments ) {
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
}

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