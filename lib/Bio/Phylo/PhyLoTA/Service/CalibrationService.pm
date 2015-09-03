package Bio::Phylo::PhyLoTA::Service::CalibrationService;
use strict;
use warnings;
use JSON;
use File::Temp 'tempfile';
use LWP::UserAgent;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Service;
use Bio::Phylo::PhyLoTA::Domain::FossilData;
use Bio::Phylo::PhyLoTA::Domain::CalibrationTable;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::Util::CONSTANT ':namespaces';
use base 'Bio::Phylo::PhyLoTA::Service';

=head1 NAME

Bio::Phylo::PhyLoTA::Service::CalibrationService - manages tree calibration

=head1 SYNOPSYS

 use Bio::Phylo::PhyLoTA::Service::CalibrationService;
 use Bio::Phylo::IO 'parse_tree';
 
 my $cs = Bio::Phylo::PhyLoTA::Service::CalibrationService->new;
 
 # read a PROTEUS compatible tab-separated spreadsheet with fossils,
 # returns an array of Bio::Phylo::PhyLoTA::Domain::FossilData objects
 my @fossils = $cs->read_fossil_table('fossils.tsv');
 
 # identify the calibration points the nodes should attach to
 my @identified = map { $cs->find_calibration_point($_) } @fossils;
 
 # read a newick tree
 my $tree = parse_tree(
 	'-format' => 'newick',
 	'-file'   => 'tree.dnd',
 	'-as_project' => 1,
 );
 
 # creates a Bio::Phylo::PhyLoTA::Domain::CalibrationTable
 my $ct = $cs->create_calibration_table( $tree, @identified );
 
 # calibrate the tree
 my $chronogram = $cs->calibrate_tree(
 	'-tree'              => $tree,
 	'-numsites'          => 1140,
 	'-calibration_table' => $ct,
 );

=head1 DESCRIPTION

This service class manages the workflow with respect to tree calibration. This includes
reading of PROTEUS-compatible fossil data sets ('read_fossil_table'), identifying the
corresponding taxa they represent ('find_calibration_points'), turning these into a
calibration table that represents the relevant calibration points given the focal tree
('create_calibration_table') and finally using these to produce a chronogram.

=head1 METHODS

=over

=item new

The constructor takes no arguments.

=cut

sub new {
    my $class = shift;
    my $self = bless {}, $class;
    return $self;
}

=item calibrate_tree

Calibrates an additive tree using a CalibrationTable object, returns an ultrametric tree.
Arguments:

 -numsites => Number of sites in the alignment that induced the tree
 -tree     => A Bio::Phylo::Forest::Tree object
 -calibration_table => A Bio::Phylo::PhyLoTA::Domain::CalibrationTable object
 -treepl_smooth     => (Optional) TreePL smoothing factor
 -nthreads          => (Optional) Number of threads

=cut

sub calibrate_tree {
	my ($self,%args) = @_;
	my $config    = $self->config;
	my $smooth    = $args{'-treepl_smooth'}     || $config->TREEPL_SMOOTH;
	my $numsites  = $args{'-numsites'}          || die "Need -numsites argument";
	my $ct        = $args{'-calibration_table'} || die "Need -calibration_table argument";
	my $tree      = $args{'-tree'}              || die "Need -tree argument";
	my $nthreads  = $args{'-nodes'}             || 1;
	
	my $seed = $config->RANDOM_SEED;
	my ( $ifh, $writetree ) = tempfile();
	my ( $ofh, $readtree )  = tempfile();
	my ( $tfh, $tplfile )   = tempfile();
	print $ofh $tree->to_newick( nodelabels => 1 );
	close $ofh;
	close $ifh;
	
	# print treePL config file header
	print $tfh <<"HEADER";
treefile = $readtree
smooth = $smooth
numsites = $numsites
outfile = $writetree
nthreads = $nthreads
seed = $seed
HEADER

	# print MRCA statements
	print $tfh $ct->to_string;
	# run treePL
	close $tfh;
	$self->logger->info("Wrote treePL config file to  $tplfile");
	system( $config->TREEPL_BIN, $tplfile ) && die $?;
	if ( -e $writetree and -s $writetree ) {        
        	my $result = parse_tree(
			'-format'     => 'newick',
			'-file'       => $writetree,
			'-as_project' => 1,
		);
		unlink $readtree, $writetree, $tplfile, "${tplfile}.r8s", "${writetree}.r8s";
        	return $result;
	}
	else {
		$self->logger->fatal("TreePL failed, reconsider your calibration points.");
#		unlink $readtree, $writetree, $tplfile;
		exit(1);
	}
}

=item find_calibration_point

Given a L<Bio::Phylo::PhyLoTA::Domain::FossilData>, this method looks at the
CalibratedTaxon from the fossil table and looks up this taxon name in the
phylota nodes table. If nodes are found, they are attached. Note that usually, 
one or no node will be found but for some taxon names there are multiple names in
the nodes table, due to homonyms across kingdoms.

=cut

sub find_calibration_point {
    my ( $self, $fd ) = @_;
    my $log = $self->logger;
    my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
    
    my @nodes;

    # search for all nodes for the calibrated taxon
	for my $ct ( @{$fd->calibrated_taxon }) {
    	push @nodes, $mts->get_nodes_for_names( $ct );
    }
    
    # if nodes are found, attach it as a calibration point 
    # to the FossilData object
    $log->debug("found ".scalar(@nodes)." nodes for fossil ".$fd->fossil_name);

    if ( scalar @nodes > 0 ){
		$fd->calibration_points(@nodes);
		for ( @nodes ) {
			$log->info("Found node " .  $_->taxon_name . " for calibrated taxon");					
		}
    }
    return $fd;
}

=item create_calibration_table

Given a tree and an array of fossils, creates a 
L<Bio::Phylo::PhyLoTA::Domain::CalibrationTable>. Note that this method assumes for now 
that whitespaces in species names in the tree are substituted to underscores 
(e.g. Mandevilla emarginata -> Mandevilla_emarginata).

For crown fossils, the mrca of the determined node will be the calibrated node.
For stem fossils, we set the node to be calibrated as the parent of the calibrated 
taxon found in the database; since there are only few extinct animals in the NCBI taxonomy, 
which are most likely not present in the inferred tree, this parent node gives the best 
approximation to the stem age. 

=cut

sub create_calibration_table {
	my ( $self, $tree, @fossildata ) = @_;
	
	my $config = $self->config;
	my $logger = $self->logger;        
	my $table  = Bio::Phylo::PhyLoTA::Domain::CalibrationTable->new;
	my $cutoff = $config->FOSSIL_BEST_PRACTICE_CUTOFF;
	$tree->set_namespaces( 'fig' => _NS_FIGTREE_ );

	# expand to all terminal taxa cf. the taxonomy
	my $treestr = $tree->to_newick;
	$logger->debug("Calibrating tree $treestr");

	my %taxa_in_tree = map { $_->get_name => 1 } @{$tree->get_terminals};
	
	my @all_mrcas; 
	FOSSIL: for my $fd ( @fossildata ) {

		# find all descendants of all calibrated higher taxa that are present in the tree
		$logger->debug("Processing fossil " . $fd->fossil_name);
		my @nodes = $fd->calibration_points; # returns higher taxa
		my $score = $fd->best_practice_score || 0;
		if ( $score < $cutoff ) {
			$logger->warn("Quality score of fossil for calibrated taxon " .  $fd->{"Calibrated_taxon"} .  " too low. Skipping.");
			next FOSSIL;
		}

		# get all descendants of the calibrated nodes that are also in our tree
		my @subtree_nodes = map { (@{$_->get_terminals}, @{$_->get_internals}) } @nodes;
		my @descendant_ids = sort map {$_->ti} @subtree_nodes; 

        # the tree seemed to get damaged sometimes by get_descendants, so below is a workaround		
        # my @descendant_ids = sort map {$_->ti} map {@{$_->get_descendants}} @nodes;		
		
		# in some cases somehow the tree gets damaged by the operation above. 
		if ( ! $tree ) { 
			$logger->warn("Tree got lost, reloading");
			$tree = parse_tree( 
				'-format' => 'newick', 
				'-string' => $treestr, 
				);
			$tree->set_namespaces( 'fig' => _NS_FIGTREE_ );
		}
		
		my @tree_nodes = map {$tree->get_by_name($_)} @descendant_ids;
		
		if ( ! scalar(@tree_nodes) ) {
			$logger->warn("Could not calibrate fossil # " . $fd->nfos . " (" . $fd->fossil_name . "). Could not find tree nodes for calibrated taxa " . join(', ', @{$fd->calibrated_taxon}) . "in tree!");
			next FOSSIL;
		}
		
		if ( scalar(@tree_nodes) == 1 ) {
			my $node = $tree_nodes[0];
			$logger->info("Fossil " . $fd->fossil_name . " has only terminal " . $node->get_name . " to calibrate. Adding sister nodes to place on crown/stem node of terminal");
			push @tree_nodes, @{ $node->get_sisters };
		}
		
		my $mrca  = $tree->get_mrca(\@tree_nodes);
		
		# for stem fossils, take the parent of the mrca
		$mrca = $mrca->get_parent if ( lc $fd->{"CrownvsStem"} eq "stem" );

		if ( ! $mrca ) {
			$logger->warn("Could not calibrate fossil # " . $fd->nfos . " (" . $fd->fossil_name . "), no mrca found!");
			next FOSSIL;
		}
		$mrca->set_meta_object( 'fig:fossil_age_min' => $fd->min_age ) if $fd->min_age;
		$mrca->set_meta_object( 'fig:fossil_age_max' => $fd->max_age ) if $fd->max_age;
		$mrca->set_meta_object( 'fig:fossil_name' => $fd->fossil_name);
		$mrca->set_meta_object( 'fig:fossil_id' => $fd->nfos);

		push @all_mrcas, $mrca;
	}
	
	# check for nested calibrated mrcas: Do not attempt to calibrate if ancestor mrca has a younger date than a possible descendant mrca
	$logger->info("Checking fossil calibration point consistancy in tree");
	for my $i ( 0..$#all_mrcas ) {
		for my $j ( 0..$#all_mrcas ) {
			my $mrca1 = $all_mrcas[$i];
			my $mrca2 = $all_mrcas[$j];
			if (  ! ($mrca1->get_meta_object( 'fig:fossil_id' ) eq $mrca2->get_meta_object ('fig:fossil_id')) ) {
				if ( $mrca1->is_ancestor_of($mrca2) ) {
					if ( my $ancestor_max =  $mrca1->get_meta_object( 'fig:fossil_age_max') and
						 my $descendant_min = $mrca2->get_meta_object( 'fig:fossil_age_min') ) {
						if ( $ancestor_max <= $descendant_min) {
							$logger->warn("Calibrated node for fossil " 
										  . $mrca1->get_meta_object( 'fig:fossil_id') 
										  . " (" 
										  . $mrca1->get_meta_object( 'fig:fossil_name') 
										  . ") is ancestor of calibrated node for fossil "
										  . $mrca2->get_meta_object( 'fig:fossil_id') 
										  . " (" 
										  . $mrca2->get_meta_object( 'fig:fossil_name') 
										  . ") although the fossil for the latter is older! Cannot create calibration table.");							
							return(0);
						}
					}			   
				}
			}
		}
	}
	
	# create row in taxa table with all leaf descendants of mrca
	for my $mrca ( @all_mrcas ) {

		my @ids = grep {/^\d+$/} map {$_->get_name} @{$mrca->get_terminals};
	
		$table->add_row(
			'taxa'    => [ sort { $a cmp $b } @ids ],
			'min_age' => $mrca->get_meta_object( 'fig:fossil_age_min'),#$fd->min_age,
			'max_age' => $mrca->get_meta_object( 'fig:fossil_age_max'),#$fd->max_age,
			'name'    => $mrca->get_meta_object( 'fig:fossil_name'),#$fd->fossil_name,	    
			'nfos'    => $mrca->get_meta_object( 'fig:fossil_id'),#$fd->nfos,	    
	        );  
	}		
	$table->sort_by_min_age;

	return $table;	
}

=item read_fossil_table

Given a fossil table, returns an array of 
L<Bio::Phylo::PhyLoTA::Domain::FossilData> objects.

=cut

sub read_fossil_table {
    my ( $self, $file ) = @_;    
    open my $fh, '<', $file or die $!;
    my ( @records, @header );
    LINE: while(<$fh>) {
        chomp;
        my @fields = split /\t/, $_;
        if ( not @header ) {
            @header = @fields;
            next LINE;
        }
        if ( @fields ) {
            my %record;
            $record{$header[$_]} = $fields[$_] for 0 .. $#header;
            push @records, Bio::Phylo::PhyLoTA::Domain::FossilData->new(\%record);
        }
    }
    # calibratedTaxon field can have more than calibrated taxon split by commata
    close $fh;

    for my $r ( @records ) {
    	my $ct = $r->{"CalibratedTaxon"};
    	my @taxa = split(',', $ct);
    	$r->{"CalibratedTaxon"} = \@taxa;
    }

    return @records;
}

=item fetch_fossil_dates

Fetches an array of L<Bio::Phylo::PhyLoTA::Domain::FossilData> from
L<http://fossilcalibrations.org>. Note that this requires a taxon name
as an argument. "Primates" appears to work. 

XXX: Note that the CalibratedTaxon field that is returned may sometimes
create difficulties for any subsequent TNRS. For example, for Primates
we get the following names:

- Primates
- Strepsirhini
- Anthropoidea
- Catarrhini
- Hominoidea
- Hominidae
- Chimpanzee-Human (!!!)
- Humanity (!!!)

=cut

sub fetch_fossil_dates {
	my ( $self, $taxon ) = @_;
	my $url = 'http://fossilcalibrations.org/api/v1/calibrations?clade=' . $taxon;
	my $ua  = LWP::UserAgent->new;
	my $res = $ua->get($url);
	if ( $res->is_success ) {
		my @records;
		my $data = decode_json $res->decoded_content;
		
		# these are the magic numbers that FC.org uses to 
		# distinguish CrownvsStem 
		my %lookup = (
			'1' => 'stem',
			'2' => 'crown',
			'3' => 'unknown',
		);
		
		# check to see if there are any results at all
		if ( ref $data->{'calibrations'} eq 'ARRAY' ) {
		
			# iterate over results
			for my $c ( @{ $data->{'calibrations'} } ) {
				delete $c->{'publicationImages'};
				delete $c->{'treeImages'};
				my $fossils = delete $c->{'fossils'};
				my $f = $fossils->[0];
				
				# here we rename keys in the hash so that they 
				# match Bio::Phylo::PhyLoTA::Domain::FossilData				
				$c->{'CalibratedTaxon'} = delete $c->{'nodeName'};
				$c->{'MinAge'} = delete $c->{'nodeMinAge'};
				$c->{'MaxAge'} = delete $c->{'nodeMaxAge'};
				$c->{'NFos'} = delete $c->{'id'};
				$c->{'CrownvsStem'} = $lookup{$f->{'locationRelativeToNode'}};
				
				# Alarm! We're parsing a literature reference.
				my $ref = delete $c->{'calibrationReference'};
				if ( $ref =~ /([^0-9]+?)\s*([0-9]{4})/ ) {
					( $c->{'NUsrcrFos'}, $c->{'DcrFos'} ) = ( $1, $2 );
				}
				push @records, Bio::Phylo::PhyLoTA::Domain::FossilData->new($c);
			}
		}
		return @records;
	}
	else {
		$self->logger->error("Couldn't fetch from fossilcalibrations.org: ".$res->status_line);
	}
}

=back

=cut

1;
