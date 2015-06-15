package Bio::Phylo::PhyLoTA::Service::CalibrationService;
use strict;
use warnings;
use File::Temp 'tempfile';
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

=cut

sub calibrate_tree {
	my ($self,%args) = @_;
	my $config    = $self->config;
	my $smooth    = $args{'-treepl_smooth'}     || $config->TREEPL_SMOOTH;
	my $numsites  = $args{'-numsites'}          || die "Need -numsites argument";
	my $ct        = $args{'-calibration_table'} || die "Need -calibration_table argument";
	my $tree      = $args{'-tree'}              || die "Need -tree argument";
	my $nthreads  = $config->NODES;
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

        $ct->remove_orphan_taxa;
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

	my %taxa_in_tree = map { $_->get_name => 1 } @{$tree->get_terminals};
	
	# find all descendants of all calibrated higher taxa that are present in the tree
	FOSSIL: for my $fd ( @fossildata ) {
		my @nodes = $fd->calibration_points; # returns higher taxa
		my $score = $fd->best_practice_score || 0;
		if ( $score < $cutoff ) {
			$logger->warn("Quality score of fossil for calibrated taxon " .  $fd->{"Calibrated_taxon"} .  " too low. Skipping.");
			next FOSSIL;
		}
			
		# expand to all terminal taxa cf. the taxonomy
		my @terminals = map { @{ $_->get_terminal_ids } } @nodes;

		# only consider terminals that are present in our tree
		@terminals = grep { $taxa_in_tree{$_} } @terminals;
		if ( ! scalar(@terminals) ) {
			$logger->warn("Could not calibrate fossil # " . $fd->nfos . " (" . $fd->fossil_name . "). Could not find tips for calibrated taxon " . $fd->calibrated_taxon . "in tree!");
			next FOSSIL;
		}

		# for stem fossils, take the parent of the mrca
		if ( lc $fd->{"CrownvsStem"} eq "stem" ) { 

			my %t = map {$_=>1} @terminals;
			my @tree_nodes =  grep { exists($t{$_->get_name})  }  @{$tree->get_terminals};					
			my $mrca = $tree->get_mrca(\@tree_nodes);
			my $parent = $mrca->get_parent;
			if ( ! $parent ) {
				$logger->warn("Could not calibrate stem fossil # " . $fd->nfos . " (" . $fd->fossil_name . ") because the mrca of all calibrated taxa has no parent node. Skipping.");
				next FOSSIL;
			}
			@terminals = map{ $_->get_name } @{$parent->get_terminals};
			@terminals =  grep { exists($taxa_in_tree{$_}) } @terminals;																							
			$mrca->set_name( $fd->fossil_name );
			$mrca->set_meta_object( 'fig:fossil_age_min' => $fd->min_age );
			$mrca->set_meta_object( 'fig:fossil_age_max' => $fd->max_age );
		}
		else {
                        # get the MRCA of the terminals, and then the tips subtended by it
                        my @tips = map { $tree->get_by_name($_) } @terminals;
                        my $mrca = $tree->get_mrca(\@tips);
                        @terminals = map { $_->get_name } @{ $mrca->get_terminals };
			$mrca->set_meta_object( 'fig:fossil_age_min' => $fd->min_age );
			$mrca->set_meta_object( 'fig:fossil_age_max' => $fd->max_age );
		}
		$table->add_row(
			'taxa'    => [ sort { $a cmp $b } @terminals ],
			'min_age' => $fd->min_age,
			'max_age' => $fd->max_age,
			'name'    => $fd->nfos,	    
	                'nfos'    => $fd->nfos,	    
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

=back

=cut

1;
