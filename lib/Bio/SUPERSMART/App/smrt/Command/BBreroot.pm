package Bio::SUPERSMART::App::smrt::Command::BBreroot;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::IO qw(parse parse_tree);

use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);


# ABSTRACT: reroots backbone tree such that the amount of paraphyletic species is minimized

=head1 NAME

BBinfer.pm - reroots backbone tree such that the amount of paraphyletic species is minimized

=head1 SYNOPSIS


=head1 DESCRIPTION

Given an input backbone tree in newick format and a list of taxa with their NCBI taxonomy identifiers
for different taxonomic ranks, re-roots the tree such that the amount of paraphyletic species
(with respect to certain taxonomic ranks) is minimized. Output is the rerooted tree in newick format.
The taxonomic level, from which the number of monophyletic subtrees is counted, is set to 
be the highest informative taxonomic level. This is, given the species list, the highest rank which 
still contains distinct taxa.

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "backbone-rerooted.dnd";
	return (
		["taxafile|t=s", "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'", { arg => "file", mandatory => 1}],
		["backbone|b=s", "a genus level backbone tree as produced by 'smrt bbinfer'", { arg => "file", mandatory => 1}],	
		["outgroup|g=s", "one or multiple taxa (names or NCBI identifiers, seperated by commata) representing the outgroup at which the tree is rerooted", {} ],
		["outfile|o=s", "name of the output tree file (in newick format), defaults to '$outfile_default'", {default => $outfile_default, arg => "filename"}],
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
		
	my $tree = parse_tree(
		'-file'   => $backbone,
		'-format' => 'newick',
	);

	# use taxon IDs instead of names
	$ts->remap_to_ti($tree);
	
	# Perform rerooting at outgroup, if given		
	if ($opt->outgroup){
		my $ogstr = $opt->outgroup;
		
		# get nodes in the tree that correspond to outgroup species and
		#  get their mrca
		my @names = split(',', $ogstr);
		# remap outgroup species names to taxon identifiers
		my @ids = map {
						(my $name = $_) =~ s/_/ /g;
						my @nodes = $ts->search_node({taxon_name=>$name})->all;
						$logger->warn("found more than one database entry for taxon $name, using first entry.") if scalar (@nodes) > 1;						
						die "could not find database entry for taxon name $name" if scalar (@nodes) == 0;						
						$nodes[0]->ti;
		} @names;
		my %og = map{ $_=>1 } @ids;
		my @ognodes = grep { exists($og{$_->get_name}) } @{$tree->get_terminals};
		my $mrca = $tree->get_mrca(\@ognodes);
		
		# reroot at mrca of outgroup
		$mrca->set_root_below;		
	} 
	
	# Try to minimize paraphyly if no outgroup given
	else {
		
		my @records = $mt->parse_taxa_file($taxafile);	
		my $level = $mt->get_highest_informative_level(@records);	
		$logger->info("highest informative taxon level : $level");
		
		# reroot the tree
		$logger->info("rerooting backone tree");	
		$tree = $ts->reroot_tree($tree, \@records, [$level]);
	}
	
	# write rerooted tree to output file
	open my $out, '>', $outfile or die $!;
	
	$ts->remap_to_name($tree);
	
	print $out $tree->resolve->to_newick;
	close $out;
	
	$logger->info("DONE, results written to $outfile");		
}

1;