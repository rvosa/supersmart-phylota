#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Tools::Run::Phylo::StarBEAST;
use Bio::Phylo::Util::Logger ':levels';

=head1 NAME

infer_clade.pl - infers species-level subtree

=head1 SYNOPSYS

 $ perl infer_clade.pl -w <directory> -n <100000> -s <10000> [--verbose] [-f <file>]

=head1 DESCRIPTION

Infers a species-level phylogeny using the multi-species, multi-locus coalescent as
implemented in *BEAST. Given a C<directory>, traverses it, looks for subdirectories
and files that match the pattern C<clade\d+/clade\d+\.xml>. These must be NeXML files
with one or more alignments that link to the same taxa block in them. Given a C<file>
argument, analyses a single file.

=cut

# process command line objects
my $verbosity = WARN;
my ( $ngens, $sfreq ) = ( 100, 100 );
my $workdir;
my $file;
GetOptions(
	'workdir=s' => \$workdir,
	'ngens=i'   => \$ngens,
	'sfreq=i'   => \$sfreq,
	'verbose+'  => \$verbosity,
	'file=s'    => \$file,
);

# instantiate helper objects
my $config = Bio::Phylo::PhyLoTA::Config->new;
my $beast  = Bio::Tools::Run::Phylo::StarBEAST->new;
my $logger = Bio::Phylo::Util::Logger->new( 
	'-level' => $verbosity, 
	'-class' => [ 
		'main',
		'Bio::Tools::Run::Phylo::StarBEAST',
		'Bio::Phylo::Parsers::Nexml',
	],
);

# configure beast
$logger->info( "setting beast executable to " . $config->BEAST_BIN );
$beast->executable( $config->BEAST_BIN );

$logger->info("setting chain length to $ngens");
$beast->chain_length($ngens);

$logger->info("setting sampling frequency to $sfreq");
$beast->sample_freq($sfreq);

$beast->beagle_SSE(1);
$beast->beagle_CPU(1);
$beast->beagle_instances(1);
$beast->overwrite(1);

# run either one file or a directory
if ( $file and -e $file ) {
	$beast->outfile_name( "${file}.beast" );
	$beast->run( $file );
}
elsif ( $workdir and -d $workdir ) {

	# iterate over entries in work dir
	opendir my $dh, $workdir or die $!;
	while( my $entry = readdir $dh ) {

		# peruse directories named cladeXXX
		if ( $entry =~ /clade\d+/ && -d "${workdir}/${entry}" ) {
	
			# this should be a nexml file with one taxa block and
			# multiple matrices
			my $file = "${workdir}/${entry}/${entry}.xml";
			if ( -e $file ) {
				$beast->outfile_name( "${file}.beast" );			
				$beast->run( $file );
			}
			else {
				$logger->warn("inconsistent directory structure, missing: $file");
			}
		}
	}
}
else {
	$logger->fatal("need -file or -workdir argument");
	exit 1;
}