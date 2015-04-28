package Bio::SUPERSMART::App::smrtutils::Command::Plot;
use strict;
use warnings;
use Template;
use List::Util 'sum';
use List::MoreUtils qw(pairwise uniq all);
use Bio::Phylo::Treedrawer;
use Bio::Phylo::IO qw(parse_tree);
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;
use Bio::Phylo::PhyLoTA::Service::CalibrationService;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);


# ABSTRACT: plots a phylogenetic tree

=head1 NAME

Plot - Visualizes a tree

=head1 SYNOPSYS

=head1 DESCRIPTION

Creates a file with a graphical representation of a tree. Currently supported formats:

=over

=item Onezoom (html format, see L<http://www.onezoom.org>)

=item jsPhyloSVG (html format, see L<http://www.jsphylosvg.org>)

=item supersmart (html format)

=back

=cut

my %formats = ( 
	'onezoom'    => 'html', 
	'supersmart' => 'html',
	'jsphylosvg' => 'html',
);

sub options {
    
	my ($self, $opt, $args) = @_;
	my $format_default    = 'supersmart';
	my $outfile_default   = "$format_default.$formats{$format_default}";
	my $width_default     = 800;
	my $height_default    = 600;
	my $style_default     = 'rectangular';
	my $bbmarkers_default = 'markers-backbone.tsv';
	my $taxa_default      = 'species.tsv';
	my $fossil_default    = 'fossils.tsv';
	my $clade_default     = 0;
	return (
		['tree|t=s', 'file name of input tree (newick format)', { arg => 'file', mandatory => 1 } ],
		['outfile|o=s', "name of the output file, defaults to '$outfile_default'", { default => $outfile_default, arg => 'file' } ],	
		['format|f=s', "output format, currently supported: OneZoom, jsPhyloSVG, SUPERSMART. Default: $format_default", { default => $format_default } ],
		['width|w=i', "visualization width. Default: $width_default", { default => $width_default } ],
		['height|h=i', "visualization height. Default: $height_default", { default => $height_default } ],
		['style|s=s', "visualization style (rectangular or circular). Default: $style_default", { default => $style_default } ],
		['markers|m=s', "backbone markers table. Default: $bbmarkers_default", { default => $bbmarkers_default } ],
		['taxa|i=s', "taxa table. Default: $taxa_default", { default => $taxa_default } ],
		['fossils|p=s', "fossil table. Default: $fossil_default", { default => $fossil_default } ],
		['clades|c=s', "Search cladeXXX folders in workdir. Default: $clade_default", { default => $clade_default } ],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		
	
	# we only have to check the 'infile' argument. 
	# if the infile is absent or empty, abort  
	my $file = $opt->tree;
	$self->usage_error('no tree argument given') if not $file;
	$self->usage_error('file $file does not exist') unless (-e $file);
	$self->usage_error('file $file is empty') unless (-s $file);			
}

sub _apply_fossil_nodes {
	my ($self,$file,$tree) = @_;
	$self->logger->info("going to apply fossil nodes using $file");
	
	# read fossil table and convert to calibration table
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $cs = Bio::Phylo::PhyLoTA::Service::CalibrationService->new;
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

sub _apply_backbone_markers {
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

sub _apply_taxon_colors {
	my ($self,$file,$tree) = @_;
	$self->logger->info("going to apply taxon colors using $file");
	my $mts = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my @records = $mts->parse_taxa_file($file);
	
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

sub _apply_clade_markers {
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

# here is my wish list for this:
# + show which taxa were exemplars (bold taxon name)
# + show the backbone branches (dashed)
# + show bootstraps as fraction, e.g. 85/100
# + show posteriors as probability 
# + show higher taxa as evolving branch colors
# + mark fossil nodes (maybe mark ranges?)
# + show branch coverage as line thickness
# + make taxa clickable (popup with link to NCBI taxonomy and to used sequences)
# - show age ranges from *BEAST
# - show bootstrapped age ranges

sub run {
	my ($self, $opt, $args) = @_;
	my $outfile = $opt->outfile;
	my $logger  = $self->logger;	
	
	# read tree, clean labels
	my $tree = parse_tree(
		'-file'   => $opt->tree,
		'-format' => 'newick',
	);
	$tree->visit(sub{
		my $n = shift;
		my $name = $n->get_name;
		$name =~ s/_/ /g;
		$name =~ s/'(.+?)'/$1/;
		$n->set_name($name);
	});
	
	# compute coordinates and round to nearest integer
	my $drawer = Bio::Phylo::Treedrawer->new(
		'-width'  => $opt->width,
		'-height' => $opt->height,
		'-tree'   => $tree,
	);
	$drawer->compute_coordinates;
	$tree->visit(sub{
		my $n = shift;
		$n->set_x( int( $n->get_x + 0.5 ) );
		$n->set_y( int( $n->get_y + 0.5 ) );
	});
	
	# apply metadata for styling
	$self->_apply_taxon_colors($opt->taxa,$tree) if -s $opt->taxa;
	$self->_apply_backbone_markers($opt->markers,$tree) if -s $opt->markers;			
	$self->_apply_fossil_nodes($opt->fossils,$tree) if -s $opt->fossils;
	$self->_apply_clade_markers($opt->clades,$tree) if -d $opt->clades;
	
	# get template file
	my $tool = lc $opt->format;
	my $file = $ENV{'SUPERSMART_HOME'} . '/data/VISUALIZATION/' . $tool . '.tmpl';
		
	# instantiate template
	my $tt = Template->new( 'ABSOLUTE' => 1 );
	my $date = localtime();
	my %args = (
		'tree'    => $tree,
		'width'   => $opt->width,
		'height'  => $opt->height,
		'style'   => $opt->style,
		'date'    => $date,
		'command' => "smrt plot @ARGV",
	);
	my $oz = $tt->process( $file, \%args, $outfile ) or die $tt->error;
	$logger->info("DONE. Tree written to $outfile");
	return 1;
}

1;
