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

 $ perl infer_clade.pl -workdir <directory> -ngens <100000> -sfreq <10000> -lfreq <10000>\
   [--verbose] [-file <file>] [--overwrite] 

=head1 DESCRIPTION

Infers a species-level phylogeny using the multi-species, multi-locus coalescent as
implemented in *BEAST. Given a C<directory>, traverses it, looks for subdirectories
and files that match the pattern C<clade\d+/clade\d+\.xml>. These must be NeXML files
with one or more alignments that link to the same taxa block in them. Given a C<file>
argument, analyses a single file.

=over

=item ngens

Specifies the number of generations. The default, which is strictly for testing, is 100.
Instead this should be something that is more on the order of 30_000_000 if not much
more.

=item sfreq

Specifies the sampling frequency. The default, which is strictly for testing, is 100.
Instead this should be something that is more on the order of 300_000 if not much
more. For example, this could be 0.1% of ngens so that you end up with a thousand trees.

=item lfreq

Species the logging frequency. The default, which is strictly for testing, is 100.
Instead this should be something that is more on the order of 300_000 if not much
more. For example, this could be 0.1% of ngens so that you end up with a thousand logged
samples.

=item overwrite

Overwrites any previously existing output files.

=cut

# process command line objects
my $verbosity = WARN;
my ( $ngens, $sfreq, $lfreq ) = ( 100, 100, 100 );
my $workdir;
my $file;
my $overwrite;
GetOptions(
	'workdir=s' => \$workdir,
	'ngens=i'   => \$ngens,
	'sfreq=i'   => \$sfreq,
	'lfreq=i'   => \$lfreq,
	'verbose+'  => \$verbosity,
	'file=s'    => \$file,
	'overwrite' => \$overwrite,
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
if ( $ngens < 30_000_000 ) {
	$logger->warn("chain length $ngens seem very low, are you just testing?");
}

$logger->info("setting sampling frequency to $sfreq");
$beast->sample_freq($sfreq);
if ( $sfreq < 300_000 ) {
	$logger->warn("sampling frequency $sfreq seem very low, are you just testing?");
}

$logger->info("setting logging frequency to $lfreq");
$beast->log_freq($lfreq);
if ( $lfreq < 300_000 ) {
	$logger->warn("logging frequency $lfreq seems very low, are you just testing?");
}

# XXX these should be configurable from phylota.ini
$beast->beagle_SSE(1);
$beast->beagle_CPU(1);
$beast->beagle_instances(1);

# overwrite any previously existing files
$beast->overwrite($overwrite);

# run either one file or a directory
if ( $file and -e $file ) {
	$beast->outfile_name( "${file}.nex" );
	$beast->logfile_name( "${file}.log" );
	$beast->run( $file );
	$logger->info("done. trees are in ${file}.nex, BEAST log in ${file}.log");
}
elsif ( $workdir and -d $workdir ) {

	# iterate over entries in work dir
	opendir my $dh, $workdir or die $!;
	while( my $entry = readdir $dh ) {

		# peruse directories named cladeXXX
		if ( $entry =~ /clade\d+/ && -d "${workdir}/${entry}" ) {
	
			# this should be a nexml file with one taxa block and
			# multiple matrices
			my $stem = "${workdir}/${entry}/${entry}";
			my $file = "${stem}.xml";
			if ( -e $file ) {
				$beast->outfile_name( "${stem}.nex" );
				$beast->logfile_name( "${stem}.log" );
				$beast->run( $file );
				my $tmpl = 'done with %s. trees are in %s.nex, log is in %s.log';
				$logger->info(sprintf $tmpl, $entry, $stem, $stem);				
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