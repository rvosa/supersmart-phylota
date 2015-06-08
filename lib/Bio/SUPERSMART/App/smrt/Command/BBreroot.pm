package Bio::SUPERSMART::App::smrt::Command::BBreroot;

use strict;
use warnings;
use List::Util 'sum';

use Bio::Phylo::Factory;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::IO qw(parse parse_tree);

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);


# ABSTRACT: reroots backbone tree

=head1 NAME

BBinfer.pm - reroots backbone tree

=head1 SYNOPSIS

 smrt bbreroot [-h ] [-v ] [-w <dir>] [-l <filename>] -t <file> -b <file> \
	[-g <names>] [-o <filename>] 

=head1 DESCRIPTION

Reroots input backbone tree(s) in newick format (--backbone or -b argument)
ideally with one or more outgroups, which can be specified as comma-separated
values for the --outgroup (or -g) argument.

If no outgroups are specified, a taxa table (as produced by smrt taxize) will
be used to attempt to reconcile the topology with the taxonomy by selecting the
rooting that minimizes the amount of taxonomic non-monophyly on the input tree.
Note that there often are multiple solutions possible and only the first one
will be (arbitrarily) selected.

A taxafile will be used regardless. By default this is the 'species.tsv' file
in the current working directory, but this can be altered using the --taxafile
or -t argument.

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "backbone-rerooted.dnd";
	my $taxa_default = "species.tsv";
	my $tree_default = "backbone.dnd";
	return (
		["taxafile|t=s", "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'", { arg => "file", default => $taxa_default}],
		["backbone|b=s", "a backbone tree(s) file as produced by 'smrt bbinfer'", { arg => "file", default => $tree_default}],	
		["outgroup|g=s", "one or multiple taxa (names or NCBI identifiers, separated by commata) representing the outgroup at "
                 ."which the tree is rerooted. Outgroup must be enclosed in quotes.", { arg => "taxon,taxon,..." } ],
		["outfile|o=s", "name of the output tree file (in newick format), defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],
		["smooth|s", "smooth tip heights left and right of root (i.e. midpointify)",{}],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	#  If taxa file or backbone tree is absent or empty, abort  
	my @files = ( $opt->taxafile, $opt->backbone );
	foreach my $file ( @files ){
		$self->usage_error("need backbone and taxafile") if not $file;
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
}

sub run {
	my ($self, $opt, $args) = @_;		
	
	# collect command-line arguments
	my $backbone = $opt->backbone;
	my $taxafile = $opt->taxafile;
	my $smooth   = $opt->smooth;
	my $outfile  = $self->outfile;
		
	# instantiate helper objects
	my $ts  = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $mt  = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $log = $self->logger;
	
	# identify outgroup once
	my $outgroup;
	if ( my $csv = $opt->outgroup ) {
		$outgroup = [map { $_->ti } $mts->get_nodes_for_names(split /,/,$csv)];
	}	
	
	# prepare taxa data
	my @records = $mt->parse_taxa_file($taxafile);
	my @ranks = ('forma', 'varietas', 'subspecies', 'species');	
	
	# open output file	
	open my $out, '>', $outfile or die $!;			
	
	# iterate over trees
	my $counter = 1;
	open my $in, '<', $backbone or die $!;
	while(<$in>) {
		my $tree = parse_tree( '-format' => 'newick', '-string' => $_ );	
	
		# use taxon IDs instead of names
		$ts->remap_to_ti($tree);
		
		# Perform rerooting at outgroup, if given		
		if ( $outgroup ) {			
			$ts->outgroup_root(
				'-tree'     => $tree,
				'-ids'      => $outgroup,
				'-ranks'    => \@ranks,
				'-records'  => \@records,
			);		
		} 
	
		# Try to minimize paraphyly if no outgroup given
		else {	
			my $level = $mt->get_highest_informative_level(@records);	
			$log->debug("highest informative taxon level : $level");		
			$tree->resolve;	
			$tree = $ts->reroot_tree($tree, \@records, [$level]);
		}
		
		# smooth the basal branch, if requested
        	if ( $smooth ) {
            		$log->debug("smoothing out diff between left and right tip heights");
			$ts->smooth_basal_split($tree);
        	}
		
		# clean up labels and write to file
		$tree = $ts->remap_to_name($tree);
		$ts->remove_internal_names($tree);
		print $out $tree->to_newick( 'nodelabels' => 1 ), "\n";
        $log->info("Rerooted backbone tree ".$counter++);		
	}
	$log->info("DONE, results written to $outfile");		
}

1;
