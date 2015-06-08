package Bio::Phylo::PhyLoTA::Service::DecorationService;
use strict;
use warnings;
use List::Util 'sum';
use List::MoreUtils qw(pairwise uniq all);
use Bio::Phylo::PhyLoTA::Service;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Service::CalibrationService;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use base 'Bio::Phylo::PhyLoTA::Service';

=head1 NAME

Bio::Phylo::PhyLoTA::Service::DecorationService - Applies decorations to trees

=head1 DESCRIPTION

Applies decorations to trees to visualize the provenance of a pipeline run, e.g.
by indicating which nodes were calibrated, to which NCBI species the tips map,
which markers were included in the supermatrix, and so on.

=over

=item apply_clade_markers

Given a directory that contains cladeXXX folders and a tree, add the clade
clusters in which each taxon participates to it, using the 'clusters' key
in the generic hash of the respective exemplar tips.

=cut


sub apply_clade_markers {
	my ($self,$dir,$tree) = @_;
	$self->logger->info("going to apply clade markers using dir $dir");
	my $mts = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my %lookup = map { $_->get_generic('species') => $_ } @{ $tree->get_terminals };

	# iterate over root directory
	opendir my $rootdh, $dir or die $!;
	while( my $entry = readdir $rootdh ) {

		# have a clade directory
		if ( $entry =~ /^clade\d+$/ and -d "${dir}/${entry}" ) {
			opendir my $cladedh, "${dir}/${entry}" or die $!;
			while( my $clade_entry = readdir $cladedh ) {

				# have a fasta file
				if ( $entry =~ /^(cluster\d+)\.fa/ ) {
					my $cluster = $1;
					my %fasta = $mts->parse_fasta_file("${dir}/${entry}/${cluster}.fa");
					my %gis = $mts->get_gis_from_fasta(%fasta);
					for my $taxon ( keys %gis ) {
						my $node = $lookup{$taxon};

						# extend list of clusters taxon participates in
						my $clusters = $node->get_generic('clusters') || [];
						push @$clusters, $cluster;
						$node->set_generic( 'clusters' => $clusters );

						# store gi(s) of focal cluster
						$node->set_generic( $cluster => $gis{$taxon} );
					}
				}
			}
		}
	}

	# collect GIs for deeper nodes
	$tree->visit_depth_first(
		'-post' => sub {
			my $n = shift;
			if ( $n->is_internal ) {

				# merge all cluster GIs from children
				my %clusters;
				for my $c ( @{ $n->get_children } ) {
					my @clusters = @{ $c->get_generic('clusters') };
					for my $cluster ( @clusters ) {
						$clusters{$cluster} = [] if not $clusters{$cluster};
						push @{ $clusters{$cluster} }, @{ $c->get_generic($cluster) };
					}
				}

				# carry forward
				%clusters = map { $_ => [ uniq @{ $clusters{$_} } ] } keys %clusters;
				$n->set_generic( 'clusters' => [ keys %clusters ], %clusters );
			}
		}
	);
}

=item apply_taxon_colors

Given a taxa table (file) and a tree, applies colors to the genera. The colors
are spaced out over the RGB spectrum from (255,127,0) to (0,127,255) (with an
oscillation in the G). The colors are blended from tips to root. To each tip,
the NCBI species ID is attached with C<set_guid>.

=cut

sub apply_taxon_colors {
	my ($self,$file,$tree) = @_;
	$self->logger->info("going to apply taxon colors using $file");
	my $mts = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my @records = $mts->parse_taxa_file($file);
	$_->{'name'} =~ s/_/ /g for @records;

	# create lookup table for genus colors, which we
	# will then blend in a traversal from tips to root
	my %genera  = map { $_->{'genus'} => 1 } @records;
	my @genera  = keys %genera;
	my $steps   = 255 / $#genera;
	for my $i ( 0 .. $#genera ) {
		my $R = int( $i * $steps );
		my $G = abs( 127 - $R );
		my $B = 255 - $R;
		$genera{$genera[$i]} = [ $R, $G, $B ];
	}

	# apply colors
	$tree->visit_depth_first(

		# pre-order
		'-pre' => sub {
			my $n = shift;
			if ( $n->is_terminal ) {

				# lookup species record and give it color of its genus
				my $name = $n->get_name;
				my ($record) = grep { $_->{'name'} eq $name } @records;
				my ( $R, $G, $B ) = @{ $genera{ $record->{'genus'} } };
				$n->set_branch_color( [ $R, $G, $B ] );
				$n->set_guid( $record->{'species'} );
				$n->set_rank( 'species' );
			}
		},

		# post-order
		'-post' => sub {
			my $n = shift;
			my @children = @{ $n->get_children };
			if ( $n->is_internal ) {

				# average over the RGB values of the children
				my @RGBs = map { $_->get_branch_color } @children;
				my $R_mean = int( sum( map { $_->[0] } @RGBs ) / @RGBs );
				my $G_mean = int( sum( map { $_->[1] } @RGBs ) / @RGBs );
				my $B_mean = int( sum( map { $_->[2] } @RGBs ) / @RGBs );
				$n->set_branch_color( [ $R_mean, $G_mean, $B_mean ] );
			}
		}
	);
}

=item apply_backbone_markers

Given a backbone marker file and a tree,
sets the following annotations:

- exemplar tips are given a true value for the 'exemplar' key in the generic hash
- exemplar tips are given a hash of marker ID to accession mappings under 'backbone'
- backbone nodes are given hash of marker ID to n seqs mappings under 'backbone'
- the root is given a hash of marker ID to marker name mappings under 'markers'

=cut

sub apply_backbone_markers {
	my ($self,$file,$tree) = @_;
	$self->logger->info("going to apply backbone markers using $file");
	my $mts = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
	my @records = $mts->parse_taxa_file($file);

	# iterate over record, add generic annotation to all nodes on the path
	# listing which markers and accessions participated in the node
	for my $r ( @records ) {
		if ( my $tip = $tree->get_by_name( $r->{'taxon'} ) ) {
			$tip->set_generic( 'exemplar' => 1 );

			# store marker ID and accession number
            my %markers;
            for my $m ( @{ $r->{'keys'} } ) {
				next if $m eq 'taxon';
                $markers{$m} = $r->{$m} if $r->{$m} and $r->{$m} =~ m/\S/;
            }
			$tip->set_generic( 'backbone' => \%markers );

			# start the counts at 1
			my %counts = map { $_ => 1 } keys %markers;
			for my $node ( @{ $tip->get_ancestors } ) {

				# may have visited node from another descendent, need
				# to add current counts
				if ( my $h = $node->get_generic('backbone') ) {
					$h->{$_} += $counts{$_} for keys %counts;
					my %copy = %$h;
					$node->set_generic( 'backbone' => \%copy );
				}
				else {
					my %copy = %counts;
					$node->set_generic( 'backbone' => \%copy );
				}
			}
		}
	}

	# get marker names
	my %markers;
	open my $fh, '<', $file or die $!;
	while(<$fh>) {
		chomp;
		if ( /# (marker\d+) cluster seed: ([^,]+),/ ) {
			my ( $marker, $acc ) = ( $1, $2 );
			my @desc = $sg->get_markers_for_accession( $acc );
			$markers{$marker} = \@desc;
			$self->logger->info("$marker ($acc): @desc");
		}
	}
	$tree->get_root->set_generic( 'markers' => \%markers );
}

=item apply_fossil_nodes

Given a fossil table and a tree, attaches a hash describing a calibration
point to each appropriate fossil node, under the 'fossil' key in the
generic hash.

=cut

sub apply_fossil_nodes {
	my ($self,$file,$tree) = @_;
	$self->logger->info("going to apply fossil nodes using $file");

	# read fossil table and convert to calibration table
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $cs  = Bio::Phylo::PhyLoTA::Service::CalibrationService->new;
	my @records = $cs->read_fossil_table($file);
	@records = map { $cs->find_calibration_point($_) } @records;

	# taxize, if need be
	my @prune;
	$tree->visit(sub{
		my $n = shift;
		my $name = $n->get_name;

		# tip doesn't have TI as name
		if ( $n->is_terminal and $name !~ /^\d+$/ ) {

			# already taxized using species.tsv
			if ( my $ti = $n->get_guid ) {
				$n->set_generic( 'binomial' => $name );
				$n->set_name($ti);
			}
			else {
				$self->logger->info("taxize: $name");
				if ( my ($node) = $mts->get_nodes_for_names($name) ) {
					$n->set_name($node->ti);
					$n->set_guid($node->ti);
					$n->set_generic( 'binomial' => $name );
				}
				else {
					push @prune, $n;
				}
			}
		}
	});
	$tree->prune_tips(\@prune) if @prune;
	my $table = $cs->create_calibration_table($tree,@records);

	# apply calibration points
	my %desc;
	$tree->visit_depth_first(
		'-post' => sub {
			my $n  = shift;
			my $id = $n->get_id;

			# start list of descendants
			if ( $n->is_terminal ) {
				$desc{$id} = [ $n->get_name ];
			}

			# extend, see if fossil found
			else {
				my @desc;
				$n->visit(sub{ push @desc, @{ $desc{shift->get_id} }});
				my @sorted = sort { $a cmp $b } @desc;
				$desc{$id} = [ uniq @sorted ];

				# iterate over fossils
				for my $f ( $table->get_rows ) {
					my @taxa = map { join ' ', split /_/, $_ } sort { $a cmp $b } $f->taxa;

					# match found
					if ( @sorted == @taxa and all { $sorted[$_] eq $taxa[$_] } 0 .. $#sorted ) {
						$n->set_generic( 'fossil' => {
							'name'    => $f->name,
							'min_age' => $f->min_age,
							'max_age' => $f->max_age,
						} );
					}
				}
			}
		}
	);
}

=back

=cut

1;
