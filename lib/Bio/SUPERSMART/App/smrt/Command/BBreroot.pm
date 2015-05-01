package Bio::SUPERSMART::App::smrt::Command::BBreroot;

use strict;
use warnings;
use List::Util 'sum';

use Bio::Phylo::Factory;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::IO qw(parse parse_tree);

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);


# ABSTRACT: reroots backbone tree

=head1 NAME

BBinfer.pm - reroots backbone tree

=head1 SYNOPSIS

smrt bbreroot [-h ] [-v ] [-w <dir>] [-l <filename>] -t <file> -b <file> [-g <names>] [-o <filename>] 

=head1 DESCRIPTION

Given input backbone tree(s) in newick format and a list of taxa with their NCBI taxonomy 
identifiers for different taxonomic ranks, re-roots the tree such that the number of 
paraphyletic taxa (with respect to certain taxonomic ranks) is minimized. Output is the 
rerooted tree(s) in newick format. The taxonomic level from which the number of 
monophyletic subtrees is counted, is set to be the highest informative taxonomic level. 
This is, given the species list, the highest rank which still contains distinct taxa.

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "backbone-rerooted.dnd";
	return (
		["taxafile|t=s", "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'", { arg => "file", mandatory => 1}],
		["backbone|b=s", "a backbone tree(s) file as produced by 'smrt bbinfer'", { arg => "file", mandatory => 1}],	
		["outgroup|g=s", "one or multiple taxa (names or NCBI identifiers, separated by commata) representing the outgroup at which the tree is rerooted. Outgroup must be enclosed in quotes.", {} ],
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
	my $outfile = $self->outfile;
	
	# instantiate helper objects
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $logger = $self->logger;
	
	# prepare taxa data
	my @records = $mt->parse_taxa_file($taxafile);
	my @ranks = ('forma', 'varietas', 'subspecies', 'species');	
	my $level = $mt->get_highest_informative_level(@records);	
	$logger->info("highest informative taxon level : $level");
	
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
		if ( $opt->outgroup ) {
		
			# parse outgroup string
			my $ogstr = $opt->outgroup;
			$logger->info("Attempting to reroot tree $counter at outgroup $ogstr");
			my @names = split(',', $ogstr);
					
			# remap outgroup species names to taxon identifiers
			my @ids = map {
				(my $name = $_) =~ s/_/ /g;
				my @nodes = $ts->search_node({ 'taxon_name' => $name })->all;
				
				# verify result
				if ( @nodes != 1 ) {
					$logger->warn("found more than one database entry for taxon $name, using first entry.") if scalar (@nodes) > 1;						
					die "could not find database entry for taxon name $name" if scalar(@nodes) == 0;
				}				
				$nodes[0]->ti;
			} @names;	
			my @all_ids = $mt->query_taxa_table(\@ids, \@ranks, @records);

			# fetch outgroup nodes and mrca
			my %og = map{ $_=>1 } @all_ids;
			my @ognodes = grep { exists($og{$_->get_name}) } @{$tree->get_terminals};
			my $mrca = $tree->get_mrca(\@ognodes);
		
			# reroot at mrca of outgroup
			$mrca->set_root_below;		
		} 
	
		# Try to minimize paraphyly if no outgroup given
		else {	
                        $tree->resolve;	
			$tree = $ts->reroot_tree($tree, \@records, [$level]);
		}
                if ( $opt->smooth ) {
                        $logger->info("smoothing out diff between left and right tip heights");
                        my $root = $tree->get_root;

                       # compute paths left and right of the root
                       my ( $left, $right ) = @{ $root->get_children };
                       my ( @hl, @hr, %h );
                       for my $tip ( @{ $left->get_terminals } ) {
                           push @hl, $tip->calc_path_to_root;
                       }
                       for my $tip ( @{ $right->get_terminals } ) {
                           push @hr, $tip->calc_path_to_root;
                       }

                       # compute averages
                       my $lm = sum(@hl)/scalar(@hl);
                       my $rm = sum(@hr)/scalar(@hr);
                       if ( $lm < $rm ) {
                           $left->set_branch_length(  $left->get_branch_length  + (abs($rm-$lm)/2) );
                           $right->set_branch_length( $right->get_branch_length - (abs($rm-$lm)/2) );
                       }
                       else {
                           $left->set_branch_length(  $left->get_branch_length  - (abs($rm-$lm)/2) );
                           $right->set_branch_length( $right->get_branch_length + (abs($rm-$lm)/2) );
                       }
                }
                $logger->info("rerooted backbone tree ".$counter++);
		
		# clean up labels and write to file
		$tree = $ts->remap_to_name($tree);
		$ts->remove_internal_names($tree);
		print $out $tree->to_newick( 'nodelabels' => 1 ), "\n";
	}
	$logger->info("DONE, results written to $outfile");		
}

1;
