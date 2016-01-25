package Bio::SUPERSMART::App::smrt::Command::Orthologize;

use strict;
use warnings;

use File::Copy qw(copy);
use File::Spec;

use Bio::SUPERSMART::Config;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;
use Bio::SUPERSMART::Service::SequenceGetter;

use Bio::SUPERSMART::Service::ParallelService;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: creates orthologous clusters of aligned sequences

=head1 NAME

Align.pm - assesses orthology in different sequence alignments and merges them into orthologous clusters

=head1 SYNOPSYS


=head1 DESCRIPTION

Given a list of aligned candidate clusters, assigns orthology among the 
clusters by performing reciprocal BLAST searches on the seed sequences 
around which the clusters were assembled. Produces a list of re-aligned 
superclusters. 

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "merged.txt";
	my $infile_default = "aligned.txt";
	return (
		[
		 "infile|i=s", 
		 "list of file locations of multiple sequence alignments  as produced by 'smrt align'", 
		 { arg => "file", default => $infile_default, galaxy_in => 1, galaxy_format => 'tabular', galaxy_type => "data" }
		],
		["outfile|o=s", 
		 "name of the output file, defaults to '$outfile_default'", 
		 {default => $outfile_default, arg => "file", galaxy_out => 1, galaxy_format => 'tabular', galaxy_type => "data" }
		],	
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	# If the infile is absent or empty, abort  
	my $file = $opt->infile;
	$self->usage_error("no infile argument given") if not $file;
	$self->usage_error("file $file does not exist") unless (-e $file);
	$self->usage_error("file $file is empty") unless (-s $file);
	
	# Outfile must be empty because we will append to it
	my $outfile = $opt->outfile;
	if ( -e $outfile and -s $outfile ) {
		$self->logger->warn("$outfile is not empty. Will overwrite.");
		open my $outfh, '>', $outfile or die $!;
		close $outfh;		
	}
}

sub run {
	my ( $self, $opt, $args ) = @_;
	
	# collect command-line arguments
	my $infile  = $opt->infile;
	my $outfile = $opt->outfile;
	my $workdir = $self->workdir;
	
	# instantiate helper objects
	my $service = Bio::SUPERSMART::Service::SequenceGetter->new;
	my $mts     = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
	my $config  = Bio::SUPERSMART::Config->new;
	my $log     = $self->logger;
	
	# parse seed GIs of aligned files
	$log->info("Going to read seed GIs from $infile");
	open my $fh, '<', $infile or die $!;
	my @gis;
	while(<$fh>) {
		chomp;
		my ( $volume,$directories,$filename ) = File::Spec->splitpath( $_ );
		if ( $filename =~ /(\d+)-.+\.fa/ ) {
			my $gi = $1;
			push @gis, $gi;
		}
	}
	close $fh;
	
	$service->merge_alignments( $config->BACKBONE_MAX_DISTANCE, $workdir, $outfile, @gis );
	
	$log->info("DONE, results written to $outfile");
	return 1;
}

1;
