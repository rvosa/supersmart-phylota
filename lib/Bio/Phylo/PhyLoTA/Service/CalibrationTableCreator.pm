# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Service::CalibrationTableCreator;
use strict;
use warnings;
use Moose;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Domain::CalibrationTable;

extends 'Bio::Phylo::PhyLoTA::Service';

=head1 NAME

Bio::Phylo::PhyLoTA::Service::CalibrationTableCreator - creates sets of calibration points

=head1 DESCRIPTION

Creates a L<Bio::Phylo::PhyLoTA::Domain::CalibrationTable> object by placing
L<Bio::Phylo::PhyLoTA::Domain::FossilData> objects on a taxonomy.

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

=item find_calibration_point

Given a L<Bio::Phylo::PhyLoTA::Domain::FossilData>, this method looks at the
CalibratedTaxon from the fossil table and looks up this taxon name in the
phylota nodes table. If nodes are found, they are attached. Note that usually, 
one or no node will be found but for some taxon names there are multiple names in
the nodes table. 

=cut

sub find_calibration_point {
    my ( $self, $fd ) = @_;
    my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
    
    #search for all nodes for the calibrated taxon
    my @nodes = $mts->get_nodes_for_names( $fd->calibrated_taxon() );
    
    #if nodes are found, attach it as a calibration point 
    #to the FossilData object
    if (scalar @nodes > 0){
	$fd->calibration_points(@nodes);
    }
    return $fd;
}



=item create_calibration_table

    Given a tree and an array of fossils, creates a 
    L<Bio::Phylo::PhyLoTA::Domain::CalibrationTable>.
    Note that this method assumes for now that 
    whitespaces in species names in the tree
    are substituted to underscores 
    (e.g. Mandevilla emarginata -> Mandevilla_emarginata).

=cut

sub create_calibration_table {
    my ( $self, $tree, @fossildata ) = @_;
    my $config = $self->config;
    my $logger = $self->logger;        
    my $table = Bio::Phylo::PhyLoTA::Domain::CalibrationTable->new;
    my $cutoff = $config->fossil_best_practice_cutoff();
    print "Cutoff : $cutoff \n";
    for my $fd ( @fossildata ) {	
	my @nodes = $fd->calibration_points();
	my $score = $fd->best_practice_score || 0;
	    if ($score >= $cutoff){	    
		if ( scalar @nodes > 0) {
		    # fetch all the tips in the taxonomy using a level-order,
		    # non-recursive, queue-based traversal
		    my @tips;
		    my @queue = @nodes;
		    while ( @queue ) {
			my $node = shift @queue;
			my @children = @{ $node->get_children };
			if ( @children ) {
			    push @queue, @children;
			}
			else {						    
			    my $name = $node->get_name;			    
			    # caution! assuming underscore-convention for 
			    # whitespaces in taxon names in the tree
			    $name =~s/\ /\_/g;	
			    # search for the equivalent node in the tree
			    if ( my $tip = $tree->get_by_name($name) ) {
				push @tips, $tip;
			    }
			}
		    }
		    if (scalar @tips > 0){		    
			# get most recent common ancestor
			my $mrca = $tree->get_mrca(\@tips);
			
			# create record in table
			$table->add_row(
			    'taxa'    => [ map { $_->get_name } @{ $mrca->get_terminals } ],
			    'min_age' => $fd->min_age,			    
			    );
		    }
		}	    
	    }
	else {
	    $logger->info("Best practice score for fossil record for taxon"
			  .$fd->calibrated_taxon()." too low.");
	}
    }
    return $table;
}

=back

=cut

1;
