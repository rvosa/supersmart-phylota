package Bio::SUPERSMART::App::smrtutils::Command::AddHigherRank;

use strict;
use warnings;

use Bio::SUPERSMART::Service::TreeService;

use Bio::Phylo::IO qw(parse);
use Bio::Phylo::Util::CONSTANT ':namespaces';

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
	my $outfile_default = 'tree-higher.nex';
	my $rank_default = 'family';
	my $format_default = 'figtree';
	return (
		['tree|t=s', 'file name of input tree (newick format)', { arg => 'file', mandatory => 1}],
		["informat|i=s", "format of input tree file, (nexus, newick, figtree) defaults to $format_default", {default => $format_default, arg => "format" }],
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
	my $rank = lc $opt->rank;
	my $logger = $self->logger;
	
	my $ts = Bio::SUPERSMART::Service::TreeService->new;

	my $tree;
	if ( lc $opt->informat  eq "figtree" ) {
		$tree = $ts->read_figtree( '-file'=> $treefile );
	} 
	else {
		$tree = parse(
			'-file'   => $treefile,
			'-format' => $opt->informat,
	    )->first;
	}
	
	my %clades;
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
			$clades{$higher_name} = [] if not $clades{$higher_name};
			push @{$clades{$higher_name}}, $n;
			$n->set_name("$higher_name-$name");
		}			
				  });
	
	# traverse clades and set name to mrca of nodes in one clade
	for my $clade ( keys %clades ) {
		my @nodes = @{$clades{$clade}};
		my $mrca = $tree->get_mrca(\@nodes);
		$logger->info("Setting internal node name $clade");
		$mrca->set_namespaces( 'fig' => _NS_FIGTREE_ );
		$mrca->set_meta_object( 'fig:name' => $clade );
	}

	$ts->write_figtree( $tree, $outfile);
	
	$logger->info("DONE. Outfile written to $outfile");
	return 1;
}

1;
