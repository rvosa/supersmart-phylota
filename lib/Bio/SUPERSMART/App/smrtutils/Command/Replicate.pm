package Bio::SUPERSMART::App::smrtutils::Command::Replicate;

use strict;
use warnings;

use File::Spec;
use Data::Dumper;

use Bio::Phylo::IO qw(parse parse_tree unparse);
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::Models::Substitution::Dna;

use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Service::ParallelService 'pfm';
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: replicate trees and alignments with simulations


=head1 NAME

Replicate - make simulated replicates of trees and alignments

=head1 SYNOPSYS

This subcommand provides functionality to generate a complete synthetic dataset from a supersmart run.
The generated synthetic dataset consists of a replicated (final) species tree and replicated alignments,
which resemble the input alignments.  
The replicated dataset can then be used as input for the SUPERSMART pipeline to validate the method and
parameters chosen for the initial run. 

=head1 DESCRIPTION

=cut

sub options {    
	my ($self, $opt, $args) = @_;
	my $aln_outfile_default = 'aligned-replicated.txt';
	my $tree_outfile_default = 'tree-replicated.dnd';
	my $taxa_outfile_default = 'taxa-replicated.tsv';
	my $format_default = 'newick';
	return (
		['tree|t=s', 'file name of input tree, must be ultrametric', { arg => 'file', mandatory => 1 } ],
		['tree_format|f=s', "format of tree input file (newick or nexus), defaults to $format_default", { default => $format_default } ],
		['alignments|a=s', "list of alignment files to replicate, as produced by 'smrt align'", { arg => 'file' } ],
		["aln_outfile|o=s", "name of output file listing the simulated alignments, defaults to '$aln_outfile_default'", {default => $aln_outfile_default, arg => "file"}],
		["tree_outfile|b=s", "name of the output tree file (newick format), defaults to '$tree_outfile_default'", {default => $tree_outfile_default, arg => "file"}],
		
		["taxa_outfile|c=s", "name the output taxa file", {default => $taxa_outfile_default, arg => "file"}],   
		["ids|i", "return NCBI identifiers in remapped tree instead of taxon names", { default=> 0}],
	    );	
}

sub validate {
	my ($self, $opt, $args) = @_;		
	
	# We only have to check the 'infile' argument. 
	#  If the infile is absent or empty, abort  
	my $file = $opt->tree;
	$self->usage_error('no tree argument given') if not $file;
	$self->usage_error("file $file does not exist") unless (-e $file);
	$self->usage_error("file $file is empty") unless (-s $file);			
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
	my $tree = parse_tree(
		'-file'   => $treefile,
		'-format' => $opt->tree_format,
	    );

	# replicate tree and write to file
	my $tree_replicated = $self->_replicate_tree($tree)->first;
	
	open my $fh, '>', $tree_outfile or die $!;	
	print $fh $tree_replicated->to_newick( nodelabels => 1 );
	close $fh;
	$logger->info("wrote replicated tree to $tree_outfile");

	$self->_write_taxafile($tree_replicated, $taxa_outfile);
	
	my @records = $mt->parse_taxa_file( $taxa_outfile );
	
	$ts->remap_to_ti( $tree, @records );
	$ts->remap_to_ti( $tree_replicated, @records );
	
	if ( my $aln = $opt->alignments ) {
		
		# read list of alignments
		$logger->info("going to read alignment file list $aln");
		open my $fh, '<', $aln or die $!;
		my @alnfiles = <$fh>;
		chomp @alnfiles; 
		close $fh;

		# some replications of alignments can take a very long time or even stall
		#  (this usually happens when replicating the binary matrix or when using phangorn's modeltest)
		#  Therefore, a timeout is set and alignments that take too long to simulate are discarded.
		my $timeout = 7200;  # set timeout to 2h
		Bio::Phylo::PhyLoTA::Service::ParallelService::timeout( $timeout );

		# replicate all alignments fiven in input align file,
		# write alignments to file and also create a list with all
		# newly written alignments
		open my $outfh, '>>', $aln_outfile or die $!; 
		
		my @replicated = pmap {			
			my ($aln) = @_; 			

			my $rep_aln = $self->_replicate_alignment( $aln, $tree, $tree_replicated );
						      		  
			if ( $rep_aln ) {
				# simulated alignment will have the same file name plus added '-simulated'
				my ( $volume, $directories, $filename ) = File::Spec->splitpath( $aln );				
				$filename =~ s/\.fa$/-replicated\.fa/g;
				$filename = $self->workdir . '/' . $filename;
				$logger->info("Writing alignment to $filename");
				unparse ( -phylo => $rep_aln, -file => $filename, -format=>'fasta' );
				print $outfh "$filename\n";
				return $filename;
			} 
			else {
				$logger->warn("Could not write replicated alignment to file; no alignment given");
			}
		} @alnfiles;
		
		close $outfh;
		$logger->info("Replicated " . scalar(@replicated) . " of " . scalar(@alnfiles) . " alignments from $aln");
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
	$logger->info("Going to make taxa table");
	my @taxa_table = $mts->make_taxa_table( \@names );
	$logger->info("Writing taxa table to $filename");
	$mts->write_taxa_file( $filename, @taxa_table );
	
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
		$logger->debug("Looking for closest non-artificial relative of artificial species $name");
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
	my ($self, $fasta, $tree, $tree_replicated) = @_; 

	my $logger = $self->logger;
	my $config = Bio::Phylo::PhyLoTA::Config->new;
	
	if ( ! ( -e $fasta and -s $fasta ) ) {
		$logger->warn("Alignment file $fasta does not exist, skipping");
		return 0;
	}
				

	$logger->info("Going to replicate alignment $fasta");
	
	my $project = parse(
		'-format'     => 'fasta',
		'-type'       => 'dna',
		'-file'     => $fasta,
		'-as_project' => 1,
	    );
	my ($matrix) = @{ $project->get_items(_MATRIX_) };

	# modeltest needs at least 3 sequences to estimate the substitution model;
	# the number of taxa in the alignment that are also present in the tree therefore
	# has to be > 2
	if ( scalar( @{$matrix->get_entities}) < 3 ) {
		$logger->warn("Cannot replicate alignment $fasta, number of taxa < 3");
		return 0;
	}
	
	$logger->info("Number of sequences in alignment $fasta : " . scalar(@{$matrix->get_entities}));

	# we have to change the definition line from something like
	# >gi|443268840|seed_gi|339036134|taxon|9534|mrca|314294/1-1140
	# to only the taxon ID: 9534
	my $matching_taxa = 0;
	my %tree_taxa = map { $_->get_name=>1 } @{ $tree_replicated->get_terminals };
	for my $seq ( @{ $matrix->get_entities } ) {
		my $seqname = $seq->get_name;
		if ( $seqname =~ m/taxon\|([0-9]+)/ ) {
			$seqname = $1;
			$logger->debug('Changing FASTA definition line from ' . $seq->get_name . " to $seqname");
			$seq->set_name($seqname);
		    }	
#		$matching_taxa++ if $tree_taxa{$seqname};
		$seq->set_generic('fasta_def_line'=>$seqname);
	}	

	# The alignment now contains as many sequences as the tree has tips.
	# We will therefore prune set of sequences. This is done by simulating a binary matrix (character is the 
	# presence/absensce of marker) using a birth-death process
	
	# collect all taxa from alignments that are present in replicated tree
	my %aln_taxa = map{ $_->get_name=>1 }  @{ $matrix->get_entities };
	my $orig_taxa_cnt =  scalar(keys(%aln_taxa));
	%aln_taxa = map {$_=>1} grep { exists $tree_taxa {$_} } keys(%aln_taxa);
		
	# it can occur that a taxon from the alignment did not end up in the final and thus also in the replicated tree.
	#  in this case, add taxa from the replicated tree to match the original number of taxa in the alignment
	while ( scalar(keys %aln_taxa) < $orig_taxa_cnt and scalar(keys %aln_taxa) < scalar(keys(%tree_taxa)) ) {
		print "Number of orig. taxa : $orig_taxa_cnt \n";
		print "Number of tree taxa : " . scalar(keys(%tree_taxa)) . "\n";
		my @terminal_ids = keys %tree_taxa;
		my $id = $terminal_ids[rand @terminal_ids];
		$logger->warn("not all taxa from alignment in replicated tree, adding random taxon $id");
		$aln_taxa{$id} = 1;
	}
	
	# make binary matrix from current alignment
	my $binary = [];
	
	for my $tax ( keys %tree_taxa ) {
		my $present = $aln_taxa{$tax} ? '1':'0';
		push @$binary, [$tax=>$present];			
	}	
	my $fac = Bio::Phylo::Factory->new;
	my $binary_matrix = $fac->create_matrix( '-matrix' => $binary);

	# make binary replicate. Here we take the relicated tree, so it is possible to have a 
	# marker for artificial species
	$logger->info("simulating binary occurence matrix for alignment $fasta");
	my $binary_rep = $binary_matrix->replicate('-tree'=>$tree_replicated, '-seed'=>$config->RANDOM_SEED);	
	if ( ! $binary_rep ) {
		$logger->warn("Cannot replicate alignment $fasta, replication of marker occurrence matrix failed");
		return 0;
	}
	$logger->debug("simulated binary matrix");

	# get taxa from replicate that have a simulated marker presence 
	my %rep_taxa =  map {$_->[0]=>1} grep {$_->[1] == 1} @{$binary_rep->get_raw};
	    	
	$logger->info(scalar(keys %aln_taxa) . " taxa in original alignment, " . scalar(keys %rep_taxa) . " taxa simulated to have marker");
	
	# it can happen that too few taxa are predicted to have the alignment,
	# in this case, randomly add taxa until the alignment will be of size 3
	if ( scalar(keys %rep_taxa) < 3 ) {
		$logger->warn("Binary simulation yielded < 3 taxa, skipping");
		return 0;
	}	

	# Lets not simulate all tips, that would take too long. Instead, 
	# prune the replicated tree to: 
	# 1. The taxa we have data for and 
	# 2. The taxa we want in our output alignment	
	my $pruned = parse('-format'=>'newick', '-string'=>$tree_replicated->to_newick)->first;
	$pruned->keep_tips( [keys %aln_taxa, keys %rep_taxa] );
	
	if ( scalar(@{$pruned->get_terminals}) < 3 ) {
		$logger->warn("Less than three taxa in pruned replicated tree. Binary simulation probably yielded too few taxa. Cannot replicate alignment $fasta.");
		return 0;
	}

	# replicate dna data: estimate model with the original tree and replicate sequences along the replicated tree
	$logger->info("Determining substitution model for alignment $fasta");
	my $model = 'Bio::Phylo::Models::Substitution::Dna'->modeltest($matrix);
	
	$logger->debug("Pruned replicated tree for sequence simulation: " . $pruned->to_newick);
	$logger->info("Simulating sequences for alignment $fasta");
	my $rep = $matrix->replicate('-tree'=>$pruned, '-seed'=>$config->RANDOM_SEED, '-model'=>$model);

	# throw out sequences that are not for our desired taxa
	for my $seq ( @{ $rep->get_entities }) {
		if ( ! $rep_taxa{$seq->get_name} ){
			$logger->debug("Removing seq for taxon " . $seq->get_name . " from replicated alignment");
			$rep->delete($seq);
			next;
		}
		# need to fix definition line to prevent a ">>" in FASTA definition line (issue #21 in Bio::Phylo)
		my $defline = $seq->get_generic('fasta_def_line');
		$defline =~ s/>//g;
		$seq->set_generic('fasta_def_line', $defline);		
	}	
	$logger->info(scalar(@{$matrix->get_entities}) . ' seqs in original, ' . scalar(@{$rep->get_entities}) . ' in replicated alignmnent');
	my @orig_names = map { $_->get_name } @{ $matrix->get_entities };
	my @rep_names = map { $_->get_name } @{ $rep->get_entities };

    return $rep;
}

1;
