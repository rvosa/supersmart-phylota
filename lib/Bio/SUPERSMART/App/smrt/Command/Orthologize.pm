package Bio::SUPERSMART::App::smrt::Command::Orthologize;

use strict;
use warnings;

use File::Copy qw(copy);
use File::Spec;

use Bio::SearchIO;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;

use Bio::Phylo::PhyLoTA::Service::ParallelService;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: creates orthologous clusters of aligned sequences

=head1 NAME

Align.pm - assesses orthology in different sequence alignments and merges them into orthologous clusters

=head1 SYNOPSYS


=head1 DESCRIPTION

Given a list of aligned candidate clusters, assigns orthology among the clusters by performing reciprocal 
blast searches on the seed sequences around which the clusters were assembled. Produces a list of 
re-aligned superclusters. 

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "merged.txt";
	my $infile_default = "aligned.txt";
	return (
		["infile|i=s", "list of file locations of multiple sequence alignments  as produced by 'smrt align'", { arg => "file", default => $infile_default}],
		["outfile|o=s", "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],	
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
	my $service = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;
	my $mts     = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $config  = Bio::Phylo::PhyLoTA::Config->new;
	my $log     = $self->logger;
	
	# parse seed GIs of aligned clusters
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
	
	# blast and cluster the seed GIs
	$log->info("Going to cluster ".scalar(@gis)." seed GIs");
	my $dbpath   = File::Spec->catfile($workdir,'seeds.fa');
	my $dbname   = $service->make_blast_db($dbpath,@gis);	
	my $report   = $service->run_blast_all($dbname);
	my @clusters = $service->cluster_blast_results($report);
	
	# merge and align
	my @cres = pmap {
		my ($clref) = @_; 
		my $id      = $clref->{'id'};
		my @gis     = @{ $clref->{'seq'} };
		my $merged  = File::Spec->catfile( $workdir, "cluster${id}.fa" );
		my $maxdist = $config->BACKBONE_MAX_DISTANCE;
					
		# turn GIs into file names 
		my @files = map { glob ( "$workdir/" .  $_ . "*.fa" ) } @gis;
		
		# profile align files to merge as many as possible
		$service->profile_align_all( $merged, $maxdist, @files );
		if ( -s $merged ) {	
			open my $outfh, '>>', $outfile or die $!;
			print $outfh $merged, "\n";
			close $outfh;
		}
	} @clusters;

	$log->info("DONE, results written to $outfile");
	return 1;
}

1;
