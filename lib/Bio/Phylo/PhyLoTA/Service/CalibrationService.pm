package Bio::Phylo::PhyLoTA::Service::CalibrationService;
use strict;
use warnings;
use File::Temp 'tempfile';
use Bio::Phylo::PhyLoTA::Service;
use Bio::Phylo::PhyLoTA::Domain::FossilData;
use Bio::Phylo::PhyLoTA::Domain::CalibrationTable;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
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

=cut

sub calibrate_tree {
	my ($self,%args) = @_;
	my $config    = $self->config;
	my $smooth    = $args{'-treepl_smooth'}     || $config->TREEPL_SMOOTH;
	my $numsites  = $args{'-numsites'}          || die "Need -numsites argument";
	my $ct        = $args{'-calibration_table'} || die "Need -calibration_table argument";
	my $tree      = $args{'-tree'}              || die "Need -tree argument";
	my $nthreads  = $config->NODES;
	my ($ifh,$writetree) = tempfile();
	my ($ofh,$readtree)  = tempfile();
	my ($tfh,$tplfile)   = tempfile();
	print $ofh $tree->to_newick;
	close $ofh;
	close $ifh;
	
	# print treePL config file header
	print $tfh <<"HEADER";
treefile = $readtree
smooth = $smooth
numsites = $numsites
outfile = $writetree
nthreads = $nthreads
HEADER

	# print MRCA statements
	print $tfh $ct->to_string, "thorough\n";

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
    
    # search for all nodes for the calibrated taxon
    my @nodes = $mts->get_nodes_for_names( $fd->calibrated_taxon() );
    
    # if nodes are found, attach it as a calibration point 
    # to the FossilData object
    $log->debug("found ".scalar(@nodes)." nodes for fossil ".$fd->fossil_name);
    if ( scalar @nodes > 0 ){
		$fd->calibration_points(@nodes);
		$log->info($_->taxon_name) for @nodes;
    }
    return $fd;
}

=item create_calibration_table

Given a tree and an array of fossils, creates a 
L<Bio::Phylo::PhyLoTA::Domain::CalibrationTable>. Note that this method assumes for now 
that whitespaces in species names in the tree are substituted to underscores 
(e.g. Mandevilla emarginata -> Mandevilla_emarginata).

=cut

sub create_calibration_table {
    my ( $self, $tree, @fossildata ) = @_;
    my $config = $self->config;
    my $logger = $self->logger;        
    my $table  = Bio::Phylo::PhyLoTA::Domain::CalibrationTable->new;
    my $cutoff = $config->FOSSIL_BEST_PRACTICE_CUTOFF;
    my $mts    = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
    
    # created mapping from ti to fossil
    $logger->info("going to map fossils to calibrated taxon IDs");
    my %fd_for_id;
    for my $fd ( @fossildata ) {
		my @nodes = $fd->calibration_points; # returns higher taxa
		my $score = $fd->best_practice_score || 0;
		if ( $score >= $cutoff ) {
			if ( scalar @nodes > 0 ) {
				for my $n ( @nodes ) {
					my $id = $n->ti;
					$fd_for_id{$id} = { 'fd' => $fd };
					
					# create mapping to descendant nodes
					my %seen;
					my @queue = @{ $n->get_children };
					while(@queue) {
						my $child = shift @queue;
						$seen{$child->ti}++;
						my @children = @{ $child->get_children };
						push @queue, @children if @children;
					}
					$fd_for_id{$id}->{'desc'} = \%seen;			
				}
			}
		} 
		else {
			$logger->info("Best practice score for fossil record for taxon"
				  .$fd->calibrated_taxon()." too low.");
		}		   	
    }
    
    # assign each tip to the most recent (i.e. least inclusive) calibrated
    # mrca that subtends it
    $logger->info("going to assign tree leaves to calibrated ancestral taxa");
    my %tips_for_id = map { $_ => [] } keys %fd_for_id;
    my @sorted = sort { scalar(keys(%{$fd_for_id{$a}->{'desc'}})) <=> scalar(keys(%{$fd_for_id{$b}->{'desc'}})) } keys %fd_for_id;
    for my $tip ( @{ $tree->get_terminals } ) {
    	my $name = $tip->get_name;
    	$logger->debug("trying to assign leaf $name...");
    	my $node;
    	if ( $name =~ /^\d+$/ ) {
    		$node = $mts->find_node($name);
    	}
    	else {
    		$name =~ s/_/ /g;
    		($node) = $mts->get_nodes_for_names($name);
    	}
    	my $id = $node->ti;
		for my $fid ( @sorted ) {
			if ( $fd_for_id{$fid}->{'desc'}->{$id} ) {
				push @{ $tips_for_id{$fid} }, $tip;
				my $ct = $fd_for_id{$fid}->{'fd'}->calibrated_taxon;
				$logger->debug("assigning leaf $name to $ct");
			}
		}
    }
    
    # add the fossils that we could place on the tree to the table
    $logger->info("going to add within-clade calibrated taxa in table");
    for my $id ( keys %tips_for_id ) {
    	if ( scalar @{ $tips_for_id{$id} } ) {
    	
			# create record in table
			$table->add_row(
				'taxa'    => [ sort { $a cmp $b } map { $_->get_name } @{ $tips_for_id{$id} } ],
				'min_age' => $fd_for_id{$id}->{'fd'}->min_age,
				'name'    => $fd_for_id{$id}->{'fd'}->calibrated_taxon,	    
			);    	
    	}
    	else {
    	
    		# this can happen, for example, if the CalibratedTaxon was a 
    		# homonym in a different kingdom
    		$logger->warn("could not find descendants for id $id in focal tree");
    	}    
    }
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
    close $fh;
    return @records;
}

=back

=cut

1;