package Bio::SUPERSMART::App::smrtutils::Command::Replicate;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse unparse);
use Bio::Phylo::Util::CONSTANT ':objecttypes';

use Bio::Phylo::PhyLoTA::Service::TreeService;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: 

=head1 NAME

Simulate - 

=head1 SYNOPSYS


=head1 DESCRIPTION

=cut

my %formats = ('onezoom' => 'htm');

sub options {
    
	my ($self, $opt, $args) = @_;
	return (
		['tree|t=s', 'file name of input tree (newick format), must be ultrametric', { arg => 'file', mandatory => 1 } ],
		['alignments|a=s', "list of alignment files to replicate, as produced by 'smrt align'", { arg => 'file' } ]
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

	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;

	# read tree
	my $tree = parse(
		'-file'   => $treefile,
		'-format' => 'newick',
	    )->first;

	# we will work with identifiers (in fasta files and trees)
	$ts->remap_to_ti( $tree );

	my $tree_replicated = $self->_replicate_tree($tree)->first;

	if ( my $aln = $opt->alignments ) {

		# read list of alignments
		$logger->info("going to read alignment file list $aln");
		open my $fh, '<', $aln or die $!;
		my @alnfiles = <$fh>;
		chomp @alnfiles; 
		close $fh;
		
                # get array of replicated alignment objects	      
		my @alns = map { $self->_replicate_alignment($_, $tree_replicated ) } @alnfiles;

		# write alignments to file
		my $aln_list = 'aligned-simulated.txt';
		open my $outfh, '>>', $aln_list or die $!;
		for my $i ( 0..$#alns ) {
			my $cnt = $i+1;
			my $filename = "simulated-aln-$cnt.fa";
			$logger->info("Writing alignment to $filename");
			unparse ( -phylo => $alns[$i], -file => $filename, -format=>'fasta' );
			print $outfh "$filename\n";
		}
		close $outfh;
	}
	
	return 1;
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
	
	$logger->info("Going to replicate alignment $fasta");
	
	my $project = parse(
		'-format'     => 'fasta',
		'-type'       => 'dna',
		'-file'     => $fasta,
		'-as_project' => 1,
	    );
	my ($matrix) = @{ $project->get_items(_MATRIX_) };
	
	# we haeve to change the definition line from something like
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

	# replicate dna data
	my $rep = $matrix->replicate($tree);				
	
	# need to fix definition line to prevent a ">>" in FASTA definition line (issue #21)
	for my $seq ( @{ $rep->get_entities }) {
		my $defline = $seq->get_generic('fasta_def_line');
		$defline =~ s/>//g;
		$seq->set_generic('fasta_def_line', $defline);		
	}	
	return $rep;
}

1;
