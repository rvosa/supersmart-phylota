package Bio::SUPERSMART::App::smrt::Command::Cladegraft;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Service::TreeService;

use Bio::Phylo::IO qw'parse parse_tree';

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: grafts the inferred clade trees on the backbone chronogram

=head1 NAME

Cladegraft.pm - Grafts the inferred clade trees on the backbone chronogram.

=head1 SYNOPSYS

smrt cladegraft [-h ] [-v ] [-w <dir>] -b <file> [-o <file>] [-c <file>]

=head1 DESCRIPTION

Combines a backbone tree of representative genera with one or more clade trees which have been
infered independently. Given a directory as argument 'workdir', traverses it, looks for subdirectories
and files that match the pattern clade*.nex. These must be NEXUS files. 
Given a single NEXUS file file as 'cladetree' argument, grafts this tree onto the backbone.
The resulting tree is exported in the NEWICK format.
    
=cut

sub options {
	my ($self, $opt, $args) = @_;
	my $outfile_default = "final.nex";
	my $tree_default = "consensus.nex";
	my $informat_default = 'nexus';
	my $outformat_default = 'nexus';
	return (                               	
  		["backbone|b=s", "backbone tree as produced by 'smrt bbinfer'", { arg => "file", default => $tree_default }],
		["outfile|o=s", "name of the output tree file (newick format) defaults to $outfile_default", { default=> $outfile_default, arg => "file"}],    	    
		["cladetree|c=s", "name of tree file for single clade (newick format) if only a single cladetree should be grafted", { arg => "file"}],
		["informat|i=s", "file format of the backbone tree (newick or nexus), defaults to $informat_default", { default => $informat_default }],
		["outformat|f=s", "file format of the final (newick or nexus), defaults to $outformat_default", { default => $outformat_default }],    	    		
    );
}

sub validate {
	my ($self, $opt, $args) = @_;		

	if ( ! $opt->backbone) {
		$self->usage_error("need backbone argument!");
	}

	if ( my  $file = $opt->backbone ){
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
	if ( my  $file = $opt->cladetree ) {
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);							
	}

	# check if there exist 'clade folders' in working directory
	if (! $opt->cladetree){
		my $workdir = $self->workdir;
		if ( ! -d glob "$workdir/clade*" ){			
			$self->usage_error("no cladetree argument given and no clade folders found in working directory $workdir");
		}
	}
	
	if ( $opt->informat !~ /^(?:newick|nexus)/i ) {
		$self->usage_error("format for the backbone can only be newick or nexus");
	}

	if ( $opt->outformat !~ /^(?:newick|nexus)/i ) {
		$self->usage_error("format for the output tree can only be newick or nexus");
	}

}

sub run{
	my ($self, $opt, $args) = @_;		
	
	# collect command-line arguments
	my $backbone = $opt->backbone;
	my $cladetree = $opt->cladetree || 0;
	my $outfile = $self->outfile;
	my $workdir = $self->workdir;
	
	# instantiate helper objects
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $logger = $self->logger;
	
	# parse backbone tree
	my $backbone_tree = parse_tree(
		'-file'   => $backbone,
		'-format' => $opt->informat,
		'-ignore_comments' => 1,
	);

	#my ( $newicktree, %map ) = $ts->parse_newick_from_nexus( $backbone, '-ignore_comments' => 0, '-id_map' => 1 );
	
	for my $n (@{$backbone_tree->get_terminals}){
		print $n->get_generic('length_median') . "\n";
		print $n->get_name . "\n";
	}

	# filenames for nexus trees written by BEAST
	$ts->remap_to_ti($backbone_tree);
	my $grafted = $backbone_tree;
	
	# graft single clade tree 
	if ( $cladetree=~ /(clade\d+)/ ) {
		$logger->info("Grafting tree $cladetree");
		my $stem = $1;
		$grafted = $self->_graft_single_tree( $grafted, $stem, $ts );
	}
	# graft all clades in working directory
	else {
	    $logger->info("Grafting all clades in working directory $workdir");
	    opendir my $dh, $workdir or die $!;
	    while( my $entry = readdir $dh ) {                        
		    if ( $entry =~ /clade\d+/ && -d "${workdir}/${entry}" ) {                                
				    $logger->info( "Processing clade $entry" );
				    $grafted = $self->_graft_single_tree( $grafted, $entry, $ts );					
		    }       
	    }
    }
    # save final tree with taxon names in newick format

    $logger->info('Retreiving taxon names for final tree');
    $ts->remap_to_name($grafted);    
    
    my $string;
    if ( $opt->outformat =~ /nexus/i ) {
	    $string = $grafted->to_nexus( 'nodelabels' => 1 );
    }
    else {
	    $string = $grafted->to_newick( 'nodelabels' => 1 );
    }
    open my $outfh, '>', $outfile or die $!;        
    print $outfh $string; 
    close $outfh;	

    $logger->info("DONE, results written to $outfile");	
    
}

sub _graft_single_tree {
	my ( $self, $tree, $clade, $ts ) = @_;

	my $logger = $self->logger;

	# make file names
	my $workdir = $self->workdir;
	# this should be a nexus file			
	my $stem = "${workdir}/${clade}/${clade}";
	my $file = "${stem}.nex";
	my $ogfile = "${workdir}/${clade}/outgroup.txt";

	if ( ! -e $file ) {
		$logger->warn("could not find file $file for clade $clade");
		return;
	}
	
	# make consensus from nexus file
	my $consensus = $ts->consense_trees( '-infile' => $file ); 

	# save consensus tree
	open my $fh, '>', $stem.".dnd" or die $!;
	print $fh $consensus->to_newick('-nodelabels' => 1);
	close $fh;

	# also save remapped consensus tree                
	my $remapped_consensus = parse('-string'=>$consensus->to_newick( nodelabels => 1 ), '-format'=>'newick')->first;
	
	$ts->remap_to_name($remapped_consensus);
	open my $fhr, '>', $stem."-remapped.dnd" or die $!;
	print $fhr $remapped_consensus->to_newick('-nodelabels' =>1);
	close $fhr;

	# prune outgroup from consensus tree (if exists)
	if ( -e $ogfile ) {
		# read outgroup from file
		my @outgroup_taxa;
		$logger->info('pruning outgroup from consensus tree');
		open my $ogfh, '<', $ogfile or die $!;
		while (<$ogfh>) {
			chomp;
			push @outgroup_taxa, $_;
		}
		close $ogfh;
		my %og = map {$_=>1} @outgroup_taxa;

		# prune outgroups from tree
		my @tipnames = map {$_->get_name} @{$consensus->get_terminals};
		@tipnames = grep {! exists($og{$_})} @tipnames;
		
		$consensus->keep_tips(\@tipnames) if scalar @tipnames;
		
		# save pruned consensus tree
		open my $fh, '>', $stem."-pruned.dnd" or die $!;
		print $fh $consensus->to_newick('-nodelabels' => 1);
		close $fh;
		
		# save pruned remapped consensus tree                
		my $remapped_pruned = parse('-string'=>$consensus->to_newick( nodelabels => 1), '-format'=>'newick')->first;
		$ts->remap_to_name($remapped_pruned);
		open my $fhr, '>', $stem."-pruned-remapped.dnd" or die $!;
		print $fhr $remapped_pruned->to_newick('-nodelabels' =>1);
		close $fhr;		
	}

	# finally graft clade tree onto backbone
	my $grafted = $ts->graft_tree( $tree, $consensus );                
		
        $grafted->visit(sub{
                my $n = shift;
                if ( $n->is_terminal ) {
			my $id = $n->get_name;
			if ($id=~m/\./){
				die("Nodename $id is decimal!");
			}
                        my $dbnode = $ts->find_node($id);
                        $logger->fatal("Could not find name for taxon id $id in database!") if not $dbnode;
			

		}
		     });
	
	return $grafted;	
}

1;
