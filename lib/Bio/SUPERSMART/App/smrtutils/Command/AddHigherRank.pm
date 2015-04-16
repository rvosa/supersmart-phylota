package Bio::SUPERSMART::App::smrtutils::Command::AddHigherRank;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Service::TreeService;

use Bio::Phylo::IO qw(parse);

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: adds information on higher taxonomic rank to all tip labels in a tree

=head1 NAME

AddHigherRank - Adds higher-level taxonomic information to the tips in a phylogenetic tree

=head1 SYNOPSYS

smrt-utils addhigherrank [-h ] [-v ] [-w <dir>] [-l <file>] [-y ] -t <file> [-o <file>] [-r ] 

=head1 DESCRIPTION

For each tip name in a given tree, looks up the classification in the NCBI taxonomy for the specified rank
(default: family) and prepends it to the tip name. 
This can be very useful to look for misplaced taxa in a tree with many tips.
Example: if a tip has name 'Coturnix coturnix' and the specified rank is 'Family',
the new tip name will be 'Phasianidae-Coturnix coturnix'.

=cut

sub options {
    
	my ($self, $opt, $args) = @_;		
	my $outfile_default = 'tree-higher.dnd';
	my $rank_default = 'family';
	return (
		['tree|t=s', 'file name of input tree (newick format)', { arg => 'file', mandatory => 1}],
		['outfile|o=s', "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => 'file'}],	
		['rank|r=s', "higher level rank to be added to tip labels, default: $rank_default", {default => $rank_default} ],		    
	    );	
}

sub validate {
	my ($self, $opt, $args) = @_;		
	
	# We only have to check the 'infile' argument. 
	#  If the infile is absent or empty, abort  
	my $file = $opt->tree;
	$self->usage_error('no tree argument given') if not $file;
	$self->usage_error('file $file does not exist') unless (-e $file);
	$self->usage_error('file $file is empty') unless (-s $file);			
}

sub run {
	my ($self, $opt, $args) = @_;    
	my $treefile = $opt->tree;
	my $outfile = $opt->outfile;
	my $rank = $opt->rank;
	my $logger = $self->logger;
	
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	
	my $tree = parse(
		'-file'   => $treefile,
		'-format' => 'newick',
	    )->first;
	
	$tree->visit( sub{
		my $n = shift;
		if ($n->is_terminal) {
			
			my $name = $n->get_name;
			$name =~ s/_/\ /g;
			$name =~ s/\'|\"//g;
			my $node = $ts->find_node({taxon_name=>$name});
			# name of the higher group, e.g. Family
			my $higher_name = '';
			
			# traverse up the tree
			my $current_rank = '';
			while ( $node and (! ( $current_rank eq $rank )) ) {
				my $tn   = $node->taxon_name;
				my $ti   = $node->ti;
				my $current_rank = $node->rank;
				$higher_name = $tn;
				$node = $node->get_parent;
				last if $current_rank eq $rank;			
			}		
			$logger->info("Changing tip name $name to $higher_name-$name");
			$n->set_name("$higher_name-$name");
		}			
		      });

	open my $out, '>', $outfile or die $!;	
	print $out $tree->to_newick(nodelabels=>1);
	close $out;
	
	$logger->info("DONE. Outfile written to $outfile");
	return 1;
}

1;
