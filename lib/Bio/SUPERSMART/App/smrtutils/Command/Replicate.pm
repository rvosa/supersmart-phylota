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

my $config = Bio::Phylo::PhyLoTA::Config->new;

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

		# replicate all alignments fiven in input align file,
		# write alignments to file and also create a list with all
		# newly written alignments
	    open my $outfh, '>', $aln_outfile or die $!;

		my @replicated = pmap {
			my ($aln) = @_;

			# set random seed to prevent issue with forking and chosing tempfile names
			my $seed = 0;
			$seed += $_ for  map {ord $_} split(//, $aln);
			srand($seed);
			
			# file name for replicated alignment
			my ( $volume, $directories, $filename ) = File::Spec->splitpath( $aln );
			$filename =~ s/\.fa$/-replicated\.fa/g;
			$filename = $self->workdir . '/' . $filename;
			$logger->debug("Checking whether replicated alignment $filename already exists");

			# replicate if not done so previously
			if ( ! -e $filename or ! -s $filename ) {
				my $rep_aln = $self->_replicate_alignment( $aln, $tree_replicated, $tree );
				
				if ( $rep_aln ) {
					# simulated alignment will have the same file name plus added '-simulated'
					$logger->info("Writing alignment $filename");
					unparse ( -phylo => $rep_aln, -file => $filename, -format=>'fasta' );
					# output average distances in alignment
					my $dist_orig = $self->_mean_dist($aln);
					my $dist_rep = $self->_mean_dist($filename);
					$logger->info("Average distance in alignments: original : $dist_orig, replicated : $dist_rep");

				}
				else {
					$logger->warn("Could not write replicated alignment to file: alignment not replicated ");
					return 0;
				}
			}
			else {
				$logger->info("Replicated alignment $filename already exists. Skipping replication.")
			}
			# write filename to alignment list
			print $outfh "$filename\n";
			
			return $filename;

		} @alnfiles;

		close $outfh;
		$logger->info("Replicated " . scalar(grep {$_} @replicated) . " of " . scalar(@alnfiles) . " alignments from $aln");
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
	my $new_ti = $mts->max_ti + 1;

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
	my ($self, $fastafile, $tree, $original_tree) = @_;

	my $logger = $self->logger;

	# create matrix object from FASTA
	$logger->fatal("Alignment file $fastafile does not exist") if  ( ! ( -e $fastafile and -s $fastafile ) );

	my $project = parse(
		'-format'     => 'fasta',
		'-type'       => 'dna',
		'-file'     => $fastafile,
		'-as_project' => 1,
	    );
	my ($matrix) = @{ $project->get_items(_MATRIX_) };
	$matrix = $self->_clean_fasta_defline( $matrix );

	if ( scalar (@{$matrix->get_entities}) < 3 ) {
		$logger->warn("Cannot replicate alignment $fastafile with less than three sequences. Skipping.");
		return 0;
	}

	# determine for which taxa we want replicated sequences
	my @rep_taxa = $self->_simulate_marker_presence( '-matrix'=>$matrix, '-tree'=>$tree, '-replace'=>0 );

	open my $fh, '>>', "all_orig_taxa.txt" or die $!;
	for my $m ( @{$matrix->get_entities} ) {
		print $fh $m->get_name . "\n";
	}
	close $fh;

	if ( scalar(@rep_taxa) > 1) {
		open my $fh, '>>', "all_rep_taxa.txt" or die $!;
		for my $t ( @rep_taxa ) {
			print $fh $t . "\n";
		}
		close $fh;
	}

	# determine substitution model for given alignment
	my $timeout = 7200; # set to 2h
	my $model = 'Bio::Phylo::Models::Substitution::Dna'->modeltest( '-matrix' => $matrix, '-timeout' => $timeout );

	# prune tree for faster sequence simulations
	my $pruned = parse('-format'=>'newick', '-string'=>$tree->to_newick)->first;
	my %keep = map {$_=>1} @rep_taxa;

	# prevent from simulating a tree with <3 tips
	my @tree_taxa = map { $_->get_name } @{ $tree->get_terminals };
	while ( scalar(keys %keep) < 3) {
		my $id = @tree_taxa[rand @tree_taxa];
		$self->logger->debug("Attempting to add taxon $id to tree taxa");
		$keep{$id} = 1;
	}
	$pruned->keep_tips( [ keys %keep ] );

	# simulate sequences
	my $rep = $matrix->replicate('-tree'=>$pruned, '-seed'=>$config->RANDOM_SEED, '-model'=>$model);

	my %rp = map {$_=>1} @rep_taxa;
	for my $r ( @{$rep->get_entities} ) {
		if ( ! $rp{$r->get_name} ) {
			$logger->info("Removing taxon " . $r->get_name . " from replicated matrix");
			$rep->delete($r);
		}
	}

	$logger->info("Number of seqs in original alignment: " . scalar(@{$matrix->get_entities}) . ", number of seqs in rep alignment: " . scalar(@{$rep->get_entities}));
	# If we had less than two simulated marker presences, the replicated alignment is not an alignment, therefore skip
	if ( @{ $rep->get_entities } < 2 ) {
		$logger->warn("Replication produced alignment with less than 2 sequences, skipping");
		return 0;
	}
	$rep = $self->_fix_fasta_defline( $rep );

	return $rep;
}

# Given an alignment as matrix object and a tree, simulates the presence in the alignment
# for each taxon in the tree using a birth-death simulation on a binary character
sub _simulate_marker_presence {
	my ($self,%args) = @_;

	my $matrix = $args{'-matrix'};
	my $tree = $args{'-tree'};

	my $logger = $self->logger;

	# collect all taxa from alignment that are also present in the tree
	my %taxa_in_aln = map{ $_->get_name=>1 }  @{ $matrix->get_entities };
	my %taxa_in_tree = map { $_->get_name=>1 } @{ $tree->get_terminals };
	my %marker_taxa = map {$_=>1} grep { exists $taxa_in_tree {$_} } keys(%taxa_in_aln);

	# if desired, taxa not present in the tree are substituted by random taxa that are in the tree
	if ( $args{'-replace'} ) {
		if ( scalar(keys %marker_taxa) < scalar(keys %taxa_in_aln) ) {
			$logger->info("Some taxa in alignment are not in the tree. Adding random taxa from tree to taxon set.");
			%marker_taxa =  $self->_add_random_taxa( \%marker_taxa, $tree, scalar(keys %taxa_in_aln) );
		}
	}

	# make binary matrix with occurrences of taxa in the alignment
	my $bin = [];
	for my $tax ( keys %taxa_in_tree ) {
		my $present = $marker_taxa{$tax} ? '1':'0';
		push @$bin, [$tax=>$present];
	}
	my $fac = Bio::Phylo::Factory->new;
	my $marker_matrix = $fac->create_matrix( '-matrix' => $bin);

	# run binary simulation to determine taxa that are simulated to be present
	my $rep_marker_matrix = $marker_matrix->replicate('-tree'=>$tree, '-seed'=>$config->RANDOM_SEED);

	# get taxa from replicate that have a simulated marker presence
	my %simulated_taxa =  map {$_->[0]=>1} grep {$_->[1] == 1} @{$rep_marker_matrix->get_raw};

	return ( keys %simulated_taxa );
}

# given a set of taxon id's and a tree, randomly add taxa from the tree to the set until it
# reaches user supplied size
sub _add_random_taxa {
	my ($self, $set, $tree, $size) = @_;

	my %taxa = %{$set};
	my @tree_taxa = map { $_->get_name } @{ $tree->get_terminals };

	while ( scalar(@tree_taxa) > scalar(keys %taxa) and scalar(keys %taxa) < $size ) {
		my $id = @tree_taxa[rand @tree_taxa];
		$self->logger->debug("attempting to add taxon $id to set of taxa");
		$taxa{$id} = 1;
	}

	return %taxa;
}

# removes excess '>' from fasta def line
sub _fix_fasta_defline {
	my ($self, $mat) = @_;
	for my $seq ( @{ $mat->get_entities }) {
		# need to fix definition line to prevent a ">>" in FASTA definition line (issue #21 in Bio::Phylo)
		my $defline = $seq->get_generic('fasta_def_line');
		$defline =~ s/>//g;
		$seq->set_generic('fasta_def_line', $defline);
	}
	return $mat;
}

# Change the FASTA definition line from the SUPERSMART style,
# e.g. ">gi|443268840|seed_gi|339036134|taxon|9534|mrca|314294/1-1140"
# to contain only the taxon ID, e.g. 9534
sub _clean_fasta_defline {
	my ($self, $matrix) = @_;

	for my $seq ( @{ $matrix->get_entities } ) {
		my $seqname = $seq->get_name;
		if ( $seqname =~ m/taxon\|([0-9]+)/ ) {
			$seqname = $1;
			$self->logger->debug('Changing FASTA definition line from ' . $seq->get_name . " to $seqname");
			$seq->set_name($seqname);
		}
		$seq->set_generic('fasta_def_line'=>$seqname);
	}

	return $matrix;
}

sub _mean_dist {
	my ($self, $filename) = @_;
	my $mt    = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	open my $fh, '<', $filename;
	read $fh, my $string, -s $fh;
	close $fh;
	my $dist = $mt->calc_mean_distance($string);
	return $dist;
}

1;
