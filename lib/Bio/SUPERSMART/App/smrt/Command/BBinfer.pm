package Bio::SUPERSMART::App::smrt::Command::BBinfer;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse parse_tree);
use Bio::Phylo::PhyLoTA::Config;
use Bio::Tools::Run::Phylo::ExaML;
use Bio::Tools::Run::Phylo::ExaBayes;
use Bio::Phylo::PhyLoTA::Service::TreeService;

use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: inference of genus-level backbone tree

=head1 NAME

BBinfer.pm - inference of genus-level backbone tree

=head1 SYNOPSIS


=head1 DESCRIPTION

Given an input supermatrix (in interleaved PHYLIP format), infers either a maximum likelihood tree 
using ExaML or uses ExaBayes to sample a posterior distribution of trees. If ExaBayes is used, a 
'consense' tree is calculated from the posterior. The tree resulting from this script is written to file.
For tree inference with examl, an NCBI classification tree (in Newick format) has to be supplied 
as a starting tree. ExaML and ExaBayes produce many intermediate checkpoint files, for which a
directory location needs to be specified. The software for tree inference (ExaML or ExaBayes) is
determined in the configuration file, but can optionally also be given as a command-line argumet.

=cut


sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "backbone.dnd";
	return (
		["supermatrix|s=s", "matrix of concatenated multiple sequece alignments as produced by 'smrt bbmerge'", { arg => "file", mandatory => 1}],	
		["starttree|t=s", "starting tree for ML tree inference as for instance produced by 'smrt classify'", { arg => "file", mandatory => 1}],
		["inferencetool|i=s", "software tool for backbone inference, defaults to ExaML", {default => 'examl', arg => "tool"}],			
		["outfile|o=s", "name of the output tree file (in newick format), defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],			

	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	#  If alignment or taxa file is absent or empty, abort  
	my @files = ( $opt->supermatrix, $opt->starttree );
	foreach my $file ( @files ){
		$self->usage_error("need supermatrix and starttree arguments") if not $file;
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
}


sub run {
	my ($self, $opt, $args) = @_;		
		
	# collect command-line arguments
	my $supermatrix = $opt->supermatrix;
	my $outfile = $opt->outfile;
	my $starttree = $opt->starttree;
	my $inferencetool = $opt->inferencetool;
	
	# instantiate helper objects
	my $config = Bio::Phylo::PhyLoTA::Config->new;
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $log = $self->logger;
	
	# adjust capitalization if necessary
	$inferencetool = lc $inferencetool;
	$inferencetool =~ s/examl/ExaML/g;
	$inferencetool =~ s/exabayes/ExaBayes/g;
	
	#my $tool = new "Bio::Tools::Run::Phylo::" . $inferencetool;
	
	
	
	
	
}

1;