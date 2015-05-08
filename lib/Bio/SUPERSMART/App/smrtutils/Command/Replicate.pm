package Bio::SUPERSMART::App::smrtutils::Command::Replicate;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse unparse);
use Bio::Phylo::Util::CONSTANT ':objecttypes';

use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Service::ParallelService;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: 

=head1 NAME

Replicate - make simulated replicates of trees and alignments

=head1 SYNOPSYS


=head1 DESCRIPTION

=cut

sub options {    
	my ($self, $opt, $args) = @_;
	my $aln_outfile_default = 'aligned-replicated.txt';
	my $tree_outfile_default = 'tree-replicated.dnd';
	my $taxa_outfile_default = 'taxa-replicated.dnd';
	return (
		['tree|t=s', 'file name of input tree (newick format), must be ultrametric', { arg => 'file', mandatory => 1 } ],
		['alignments|a=s', "list of alignment files to replicate, as produced by 'smrt align'", { arg => 'file' } ],
		["aln_outfile|o=s", "name of output file listing the simulated alignments, defaults to '$aln_outfile_default'", {default => $aln_outfile_default, arg => "file"}],
		["tree_outfile|f=s", "name of the output tree file (newick format), defaults to '$tree_outfile_default'", {default => $tree_outfile_default, arg => "file"}],   
		["taxa_outfile|f=s", "name the output taxa file", {default => $taxa_outfile_default, arg => "file"}],   
		["ids|i", "return NCBI identifiers in remapped tree instead of taxon names", { default=> 0}],
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

	my $logger = $self->logger;
	my $treefile = $opt->tree;
	my $aln_outfile = $opt->aln_outfile;
	my $tree_outfile = $opt->tree_outfile;
	my $taxa_outfile = $opt->taxa_outfile;
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $mt  = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;

	# read tree
	my $tree = parse(
		'-file'   => $treefile,
		'-format' => 'newick',
	    )->first;

	# replicate tree and write to file
	my $tree_replicated = $self->_replicate_tree($tree)->first;

	open my $fh, '>', $tree_outfile or die $!;	
	print $fh $tree_replicated->to_newick( nodelabels => 1 );
	close $fh;
	$logger->info("wrote replicated tree to $tree_outfile");

	$logger->info('going to write taxa file');
	$self->_write_taxafile($tree_replicated, $taxa_outfile);
	
	my @records = $mt->parse_taxa_file( $taxa_outfile );
	
	$ts->remap_to_ti( $tree_replicated, @records );
	
	if ( my $aln = $opt->alignments ) {
		
		# read list of alignments
		$logger->info("going to read alignment file list $aln");
		open my $fh, '<', $aln or die $!;
		my @alnfiles = <$fh>;
		chomp @alnfiles; 
		close $fh;

		# replicate all alignments fiven in input align file,
		# write alignments to file and also create a list with all
		# newly written alignments
		open my $outfh, '>>', $aln_outfile or die $!;		
		pmap {
			my ($aln) = @_; 
			my $rep_aln = $self->_replicate_alignment($aln, $tree_replicated );
			
			if ( $rep_aln ) {
				# simulated alignment will have the same file name plus added '-simulated'
				my $filename = $aln;
				$filename =~ s/\.fa$/-replicated\.fa/g;			
				$logger->info("Writing alignment to $filename");
				unparse ( -phylo => $rep_aln, -file => $filename, -format=>'fasta' );
				print $outfh "$filename\n";
			} 
			else {
				$logger->warn("Could not write replicated alignment to file; no alignment given");
			}
			
		} @alnfiles;
		
		close $outfh;
	}
	$logger->info("DONE. Tree written to $tree_outfile, alignment list written to $aln_outfile, taxa table written to $taxa_outfile" );
	return 1;
}

# This produces a taxa file givan a tree
sub _write_taxafile {	
	my ($self, $tree, $filename) = @_; 
	
	my $logger = $self->logger;
	my $mts    = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $mt  = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	
	my @names = map { s/_/ /g; $_;} map {$_->get_name} @{$tree->get_terminals};	
	
	# identify terminals with artificial names and 
	my %artificial = map {$_=>1} grep {! $mts->find_node({taxon_name=>$_})} @names; 	
	@names = grep {! $artificial{$_}} @names;
		
	# write initial taxa table to file
	$mts->write_taxa_file( $filename, @names );	
	
	# for the species with artificial taxon names, do as follows: first: assign an artificial
	# ti to the species. Then iteratively select the children of the father of the focal node, and
	# check if any of its children have a classification, until we found a (maybe far away) sister 
	# node that is classified. Then copy the genus, family etc. and append the entries to the table	
	my $new_ti = $mts->max_ti;
	
	my @levels = reverse $mts->get_taxonomic_ranks;

	my @records = $mt->parse_taxa_file($filename);
	
	open my $fh, '>>', $filename or die $!;
	# iterate over all artificial species, assign taxon id and append entry to taxa table
	for my $an ( keys %artificial ) {				
		(my $name = $an) =~ s/ /_/g; 
		$logger->debug("Looking for closest non-artificial relative of atrificial species $name");
		my ($node)  = grep {$_->get_name eq $name} @{$tree->get_terminals};
		my $relative;
		
		# search for the closest 'real' relative
		while ( ! $relative and ! $node->is_root ) {
			my @relatives;
			$node->visit_depth_first(-in => sub {
				my $n = shift;
				(my $current_name = $n->get_name) =~ s/_/ /g;
				if ( $current_name and (! $artificial{$current_name} )) {
					push @relatives, $current_name;					
				}
				
						 });
			$relative = $relatives[0] if scalar(@relatives);
			# if the whole clade consisted of artificial species, go one up
			$node = $node->get_parent;
		}		
		$logger->debug("Closest non-artificial relative: $relative");						
		
		# get all taxon IDs for the closest nonartificial relative, sorted from low to high ranks
		my $tid = $mts->find_node({taxon_name=>$relative})->ti;
		my ($all_ids) = grep { my %h = %{$_}; grep{ $tid } values%h } @records;

		# assign taxon ID to artificial species and take all other taxon ids from closest real relative
		$all_ids->{'species'} = $new_ti++;
		my @ids = @$all_ids{@levels};
		$logger->debug("Assinging taxon id $new_ti to artificial species $name");

		# append information to taxon table
		$name =~ s/_/ /g;
		print $fh $name . "\t" . join("\t", @ids) . "\n";
	}
	close $fh;
	$logger->info("Added " . scalar(keys %artificial) . " artificial species to taxa file $filename") if keys %artificial;
}

sub _replicate_tree {
	my ($self, $tree) = @_; 
		
	$tree->generize( 
		'-delim'     => '_', 
		'-monotypic' => 1,
		'-polypara'  => 1,
	    );
	
	my $rep = $tree->replicate( '-genera' => 1);

	return $rep;	
}

sub _replicate_alignment {
	my ($self, $fasta, $tree) = @_; 
	
	my $logger = $self->logger;
	my $config = Bio::Phylo::PhyLoTA::Config->new;
	
	$logger->info("Going to replicate alignment $fasta");
	
	my $project = parse(
		'-format'     => 'fasta',
		'-type'       => 'dna',
		'-file'     => $fasta,
		'-as_project' => 1,
	    );
	my ($matrix) = @{ $project->get_items(_MATRIX_) };
	
	$logger->info('Number of sequences in alignment $fasta : ' . scalar(@{$matrix->get_entities}));

	# we have to change the definition line from something like
	# >gi|443268840|seed_gi|339036134|taxon|9534|mrca|314294/1-1140
	# to only the taxon ID: 9534
	my $matching_taxa = 0;
	my %tree_taxa = map { $_->get_name=>1 } @{ $tree->get_terminals };
	for my $seq ( @{ $matrix->get_entities } ) {
		my $seqname = $seq->get_name;
		if ( $seqname =~ m/taxon\|([0-9]+)/ ) {
			$seqname = $1;
			$logger->debug('Changing FASTA definition line from ' . $seq->get_name . " to $seqname");
			$seq->set_name($seqname);
		    }	
		$matching_taxa++ if $tree_taxa{$seqname};
		$seq->set_generic('fasta_def_line'=>$seqname);
	}	
	# modeltest needs a tree with at least 3 species to estimate the substitution model;
	# the number of taxa in the alignment that are also present in the tree therefore
	# has to be > 2
	if ( $matching_taxa < 3 ) {
		$logger->warn("Cannot replicate alignment $fasta, number of taxa < 3");
		return 0;
	}


	# The alignment now contains as many sequences as the tree has tips.
	# We will therefore prune set of sequences. This is done by simulating a binary matrix (character is the 
	# presence/absensce of marker) using a birth-death process
	
	# make binary matrix from current alignment
	my $binary = [];
	my %aln_taxa = map{ $_->get_name=>1 } @{ $matrix->get_entities };
	for my $tax ( keys %tree_taxa ) {
		my $present = $aln_taxa{$tax} ? '1':'0';
		push $binary, [$tax=>$present];			
	}	
	my $fac = Bio::Phylo::Factory->new;
	my $binary_matrix = $fac->create_matrix( '-matrix' => $binary);

	# make binary replicate
	my $binary_rep = $binary_matrix->replicate($tree, $config->RANDOM_SEED);	
	
	# get taxa from replicate that have a simulated marker presence 
	my %rep_taxa =  map {$_->[0]=>1} grep {$_->[1] == 1} @{$binary_rep->get_raw};
	
	# it can happen that no taxon is predicted to have the alignment,
	# in this case, take the original taxa from the alignment
	if ( ! keys(%rep_taxa) ) {
		$logger->warn("setting taxa for replicating $fasta to original alignment taxa, no other taxa were predicted in binary simulation");
		%rep_taxa = %aln_taxa;
	}
	
	# Lets not simulate all tips, that would take too long. Instead, 
	# prune the tree to: 
	# 1. The taxa we have data for and 
	# 2. The taxa we want in our output alignment	
	my $pruned = parse('-format'=>'newick', '-string'=>$tree->to_newick)->first;
	$pruned->keep_tips( [keys %aln_taxa, keys %rep_taxa] );
	
	# replicate dna data
	my $rep = $matrix->replicate($pruned, $config->RANDOM_SEED);					
	# throw out sequences that are not for our desired taxa
	for my $seq ( @{ $rep->get_entities }) {
		if ( ! $rep_taxa{$seq->get_name} ){
			$rep->delete($seq);
			next;
		}
		# need to fix definition line to prevent a ">>" in FASTA definition line (issue #21 in Bio::Phylo)
		my $defline = $seq->get_generic('fasta_def_line');
		$defline =~ s/>//g;
		$seq->set_generic('fasta_def_line', $defline);		
	}	
       	$logger->info('Number of sequences in replicated alignment : ' . scalar(@{$rep->get_entities}));
	
	return $rep;
}

1;
