#!/usr/bin/perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Bio::Phylo::IO 'parse';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Tools::Run::Phylo::ExaML;

=head1 NAME

infer_backbone.pl - infers backbone tree using ExaML

=head1 SYNOPSIS

 $ infer_backbone.pl -o <file> -w <dir> -c <common tree> -s <supermatrix> [--verbose]

=head1 DESCRIPTION

Given an input supermatrix (in interleaved PHYLIP format) and an NCBI classification
tree (in Newick format), infers a maximum likelihood tree using ExaML and writes this
to an output file. ExaML produces many intermediate checkpoint files, for which a
directory location needs to be specified.

=cut

# process command line arguments
my $verbosity = WARN;
my ( $outfile, $workdir, $commontree, $supermatrix );
GetOptions(
	'outfile=s'     => \$outfile,
	'workdir=s'     => \$workdir,
	'commontree=s'  => \$commontree,
	'supermatrix=s' => \$supermatrix,
	'verbose+'      => \$verbosity,
);

# instantiate helper objects
my $config = Bio::Phylo::PhyLoTA::Config->new;
my $examl  = Bio::Tools::Run::Phylo::ExaML->new;
my $logger = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# set outfile name
$logger->info("going to create output file $outfile");
$examl->outfile_name($outfile);

# set working directory
$logger->info("going to use working directory $workdir");
$examl->work_dir($workdir);

# set examl location
$logger->info("going to use ExaML executable ".$config->EXAML_BIN);
$examl->executable($config->EXAML_BIN);

# set parser location
$logger->info("going to use parser executable ".$config->PARSER_BIN);
$examl->parser($config->PARSER_BIN);

# set mpirun location
$logger->info("going to use mpirun executable ".$config->MPIRUN_BIN);
$examl->mpirun($config->MPIRUN_BIN);

# set number of nodes
$logger->info("setting number of MPI nodes ".$config->NODES);
$examl->nodes($config->NODES);

# set substitution model
$logger->info("setting substitution model ".$config->EXAML_MODEL);
$examl->m($config->EXAML_MODEL);

# set examl to quiet
$logger->info("setting ExaML to requested verbosity");
$examl->quiet($verbosity <= WARN);

# here we need to read the names from the phylip file and adjust the
# common tree accordingly: it must only retain the taxa in the supermatrix,
# and it must be fully resolved.
my $intree = File::Spec->catfile( $workdir, 'user.dnd' );
open my $fh, '>', $intree or die $!;
print $fh parse(
	'-format' => 'newick',
	'-file'   => $commontree,
)->first
	->keep_tips( [ read_tipnames($supermatrix) ] )
	->resolve
	->remove_unbranched_internals
	->deroot
	->to_newick;

# run the analysis
$logger->info("patience please, running $supermatrix with starting tree $intree");
my $backbone = $examl->run(
	'-phylip' => $supermatrix,
	'-intree' => $intree,
);
$logger->info("done, result written to $backbone");

# reads the supermatrix (phylip format) and returns the tip names from it
sub read_tipnames {
	my $supermatrix = shift;
	my $ntax = 0;
	my $line = 0;
	my @result;	
	open my $fh, '<', $supermatrix or die $!;
	LINE: while(<$fh>) {
		chomp;
		my $word;
		if ( /^(\S+)/ ) {
			$word = $1;
		}
		if ( not $ntax ) {
			$ntax = $word;
			$logger->debug("$supermatrix has $ntax taxa");
			next LINE;
		}
		push @result, $word;
		$logger->debug("adding taxon $word");
		last LINE if ++$line == $ntax;
	}
	return @result;
}
