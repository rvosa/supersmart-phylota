package Bio::SUPERSMART::App::smrtutils::Command::SuggestOutgroup;

use strict;
use warnings;

use List::MoreUtils qw(uniq);

use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);



# ABSTRACT: suggests an outgroup for a group of species

=head1 NAME

SuggestOutgroup.pm - Lists possible outgroups for a list of taxa or 

=head1 SYNOPSYS

smrt-utils suggestoutgroup [-h ] [-v ] [-w <dir>] [-l <file>] [-y ] [-n <file>] [-n <file>] [-r ] [-g ] 

=head1 DESCRIPTION

=cut

sub options {
    
	my ($self, $opt, $args) = @_;		
	my $rank_default = 'genus';
	return (
		[ 'infile|n=s',  'file containing a list of taxon names', { arg => 'file'}],
		[ 'taxafile|t=s',  "taxa file as produced by 'smrt taxize'", { arg => 'file'}],
		[ 'root_taxon|r=s', "higher-level taxon",{} ],		
		[ 'outgroup_rank|g=s', "rank of outgroups, defaults to $rank_default", { default=>$rank_default } ],		
		
	    );	
}

sub validate {
	my ($self, $opt, $args) = @_;		
	$self->usage_error('need either --names, --taxafile  or --root_taxon argument ') if not ( $opt->infile || $opt->root_taxon || $opt->taxafile);
}

sub run {
	my ($self, $opt, $args) = @_;    
	
	my $log = $self->logger;
	my $mts = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;

	my $outgroup_rank = lc $opt->outgroup_rank;

	my $root_node;

	# parse input taxa and get nodes
	if ( my $infile = $opt->infile ) {
		my @names;
		open my $fh, '<', $infile or die $!;
		@names = <$fh>;
		chomp(@names);
		# remove empty lines if present
		@names = grep /\S/, @names;
		close $fh;
		$log->info( "read " . scalar(@names) . " species names from $infile" );
		
		my @records = $mts->make_taxa_table(@names);

		my $root_rank = $mt->get_root_taxon_level( @records );
		my $root_ti = $records[0]{$root_rank};
		
		$log->debug("Root taxon id : $root_ti");
		
		# get node from database and search for siblings and their descendants
		$root_node = $mts->find_node({ti=>$root_ti});
	}
	
	elsif ( my $root_taxon = $opt->root_taxon ) {
		$log->debug("Root taxon $root_taxon given as input");
		$root_node = $mts->search_node({taxon_name=>$root_taxon})->first;
	}
	elsif ( my $taxafile = $opt->taxafile ) {
		my @records = $mt->parse_taxa_file($taxafile);
    
		my $root_rank = $mt->get_root_taxon_level( @records );
		
		my $root_ti = $records[0]{$root_rank};		
		$log->debug("Root taxon id : $root_rank");		
	}
		
	# get siblings of the root node and determine all descendents of rank of interest
	my @siblings = @{$root_node->get_siblings};
	my @sibling_descendants = map{ @{ $_->get_descendants_at_rank($outgroup_rank)}  }  @{$root_node->get_siblings};
	
	my @results;
	for my $og_node ( @sibling_descendants ) {

		$log->debug( "Possible outgroup  : " . $og_node->taxon_name . "( $outgroup_rank )");
		
		# get all terminals for outgroup
		my @terminal_ids = @{$og_node->get_terminal_ids};
		next if not scalar (@terminal_ids);

		my @terminal_nodes = map { $mts->find_node({ti=>$_}) } @terminal_ids;
		my @ranks = uniq ( map{ $_->rank } @terminal_nodes) ;

		# get number of sequences for whole outgroup
		my @seq_cnt = map { scalar($mts->search_seq({ti=>$_})) } @terminal_ids; 
		
		my $seqs_for_og;
		$seqs_for_og += $_ for @seq_cnt;

		my $str = "\tMembers (" . join (',', @ranks ) . ") : " . scalar(@terminal_ids) . " total seqs :  $seqs_for_og ";
		$log->debug($str);
		
		my %res;
		$res{'name'} = $og_node->taxon_name;
		$res{'members'} = scalar(@terminal_ids);
		$res{'types'} = \@ranks;
		$res{'seqs'} = $seqs_for_og;
		$res{'coverage'} = $seqs_for_og / scalar(@terminal_ids);
		push @results, \%res;
	}
	
	# sort by coverage and remove monotypic groups
	@results = grep { $_->{'members'} > 1 } @results;
	my @sorted_results = sort {  $b->{'coverage'} <=> $a->{'coverage'} } @results;
	
	
	$log->info("Found " . scalar(@sorted_results) . " possible outgroups (sorted by sequence coverage): ");
	for my $r ( @sorted_results ) {
		$log->info("$outgroup_rank:\t" . $r->{'name'} . ",\tmembers: " . $r->{'members'} . ",\taverage number of sequences per member: " . sprintf "\t%.2f", $r->{'coverage'});
	}
}

1;
