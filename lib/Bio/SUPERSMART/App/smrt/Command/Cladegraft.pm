package Bio::SUPERSMART::App::smrt::Command::Cladegraft;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Util::CONSTANT qw':objecttypes :namespaces';
use Bio::Phylo::IO qw'parse parse_tree';

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

my $conf = Bio::Phylo::PhyLoTA::Config->new;

# ABSTRACT: grafts the inferred clade trees on the backbone chronogram

=head1 NAME

Cladegraft.pm - Grafts the inferred clade trees on the backbone chronogram.

=head1 SYNOPSYS

smrt cladegraft [-h ] [-v ] [-w <dir>] -b <file> [-o <file>] [-c <file>]

=head1 DESCRIPTION

Combines a backbone tree of representative genera with one or more clade trees 
that have been inferred independently. Given a directory as argument 'workdir', 
traverses it, looks for subdirectories and files that match the pattern clade*.nex. 
These must be NEXUS files. Given a single NEXUS file file as 'cladetree' argument, 
grafts this tree onto the backbone. The resulting tree is exported in the NEXUS 
(figtree) format.
    
=cut

sub options {
	my ($self, $opt, $args) = @_;
	my $outfile_default = "final.nex";
	my $tree_default    = "consensus.nex";
	my $heights_default = $conf->NODE_HEIGHTS;
	my $squish_default  = 'none';
	return (                               	
  		["backbone|b=s", "backbone tree as produced by 'smrt consense'", { arg => "file", default => $tree_default }],
		["outfile|o=s", "name of the output tree file (newick format) defaults to $outfile_default", { default=> $outfile_default, arg => "file"}],    	    
		["cladetree|c=s", "name of tree file for single clade (newick format) if only a single cladetree should be grafted", { arg => "file"}],
		["heights|e=s", "node heights (ca, keep, median, mean)", { default => $heights_default, arg => 'keep|median|mean|ca' } ],
		["squish|s=s", "how to treat negative branches (zero, yulish, none)", { default => $squish_default, arg => 'zero|yulish|none' } ],
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

	if ( $opt->heights !~ /^(?:ca|keep|median|mean)$/ ) {
		$self->usage_error("heights must be one of ca, keep, median or mean");
	}
	
	if ( $opt->squish !~ /^(?:zero|none|yulish)$/ ) {
		$self->usage_error("squish must be one of yulish, zero or none");
	}
}

sub run {
	my ($self, $opt, $args) = @_;		
	
	# collect command-line arguments
	my $backbone  = $opt->backbone;
	my $cladetree = $opt->cladetree || 0;
	my $outfile   = $self->outfile;
	my $workdir   = $self->workdir;
	
	# instantiate helper objects
	my $ts      = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $logger  = $self->logger;	
	my $grafted = $ts->read_figtree( '-file' => $backbone );
	$logger->info("Read backbone $backbone");
	
	# graft single clade tree 
	if ( $cladetree=~ /(clade\d+)/ ) {
		$logger->info("Grafting tree $cladetree");
		my $stem = $1;
		$grafted = $self->_graft_single_tree( $grafted, $stem, $ts, $opt );
	}
	
	# graft all clades in working directory
	else {
	    $logger->info("Grafting all clades in working directory $workdir");
	    opendir my $dh, $workdir or die $!;
	    while( my $entry = readdir $dh ) {                        
		    if ( $entry =~ /clade\d+/ && -d "${workdir}/${entry}" ) {                                
				$logger->info( "Processing $entry" );
				$grafted = $self->_graft_single_tree( $grafted, $entry, $ts, $opt );					
		    }       
	    }
	}
	
	# map taxon IDs back to names and write to file
	$ts->write_figtree( $grafted, $outfile );
	$logger->info("DONE, results written to $outfile");	    
}

sub _graft_single_tree {
	my ( $self, $tree, $clade, $ts, $opt ) = @_;
	my $logger = $self->logger;

	# make file names
	my $workdir = $self->workdir;
	my $file    = "${workdir}/${clade}/${clade}.nex";
	my $ogfile  = "${workdir}/${clade}/outgroup.txt";
	if ( ! -e $file ) {
		$logger->warn("could not find file $file for clade $clade");
		return;
	}
	
	# make consensus from BEAST clade output
	my $consensus = $ts->consense_trees( 
		'-infile'  => $file,
		'-heights' => $opt->heights,
	);
	$consensus->set_namespaces( "fig" => _NS_FIGTREE_ );
	$consensus->get_root->set_meta_object( "fig:clade" => $clade );

	print $consensus->to_newick;
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
		$consensus->prune_tips(\@outgroup_taxa);
	}
	
	# write remapped consensus in clade directory
	my $clone = $ts->read_figtree( '-string' => $ts->to_figtree($consensus));
	my $fname    = "${workdir}/${clade}/${clade}-consensus.dnd";
 	$ts->write_figtree( $clone, $fname );

	return $ts->graft_tree( $tree, $consensus, $opt->squish );
}

1;
