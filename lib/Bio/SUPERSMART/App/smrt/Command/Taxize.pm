package Bio::SUPERSMART::App::smrt::Command::Taxize;

use strict;
use warnings;

use Bio::SUPERSMART::Config;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: writes taxa table for given list of taxon names or root taxa

=head1 NAME

Taxize.pm - reconciles taxon names with information from NCBI taxonomy

=head1 SYNOPSYS

$ smrt taxize [-v ] [-w <dir>] -i <file> [-o <file>] [-e <rank>] 

=head1 DESCRIPTION

Performs taxonomic name resolution by mapping a provided list of input names onto the 
NCBI taxonomy. Optionally expands the higher input taxa to the specified lower taxon, 
e.g. to expand a named Order to its constituent species. Produces a table that lists 
all resolved (and expanded) taxa and their higher classification. Optionally, non-binomial 
species names can be filtered out. By default, taxize also filters out taxon names containing
the keywords 'unidentified' and 'environmental_sample'.


=cut

sub options {
	my ( $self, $opt, $args ) = @_;
	my $outfile_default = "species.tsv";
	return (
		[
			"infile|i=s",
			"file with list of taxon names",
			{ arg => "file", galaxy_in => 1, galaxy_type => "file" }
		],
		[
			"outfile|o=s",
			"name of the output file, defaults to '$outfile_default'",
			{ default => $outfile_default, arg => "file", galaxy_out => 1, galaxy_type => "file" }
		],
		[
			"root_taxa|r=s",
			"one or multiple taxon names (seperated by commata) to be expanded to either 'Species' or to taxonomic rank "
			. "given in <expand_rank>. Taxon names containing spaces must be enclosed in quotes",
			{ arg => "taxon,taxon,...", galaxy_in => 1, galaxy_type => "text" }
		],
		[
			"expand_rank|e=s",
			"rank to which root taxa are expanded",
			{ default => 0, arg => "rank", galaxy_in => 1, galaxy_type => "text" }
		],
		[
			"binomials_only|b",
			"rank to which root taxa are expanded",
			{ default => 0, galaxy_in => 1, galaxy_type => "boolean" }
		],
	);
}

sub validate {
	my ( $self, $opt, $args ) = @_;

	# We only have to check the 'infile' argument.
	#  If the infile is absent or empty, abort
	my $file = $opt->infile;
	my $root = $opt->root_taxa;
	$self->usage_error("need infile and/or root_taxon argument") if not $file and not $root;
	if ( $file ) {
		$self->usage_error("file $file does not exist") unless ( -e $file );
		$self->usage_error("file $file is empty")       unless ( -s $file );
	}
}

sub run {
	my ( $self, $opt, $args ) = @_;

	# collect command-line arguments
	my $infile      = $opt->infile;
	my $expand_rank = $opt->expand_rank;
	my $root_taxa   = $opt->root_taxa;

	# instantiate helper objects
	my $log    = $self->logger;
	my $config = Bio::SUPERSMART::Config->new;
	my $mts    = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;

	my @names;

	if ( $infile ) {
		# read names from file or STDIN, clean line breaks
		open my $fh, '<', $infile or die $!;
		while(<$fh>) {

			# strip line breaks and leading/trailing whitespace
			chomp;
			s/^\s*//;
			s/\s*$//;
			push @names, $_ if /\S/;
		}
		close $fh;
		$log->info( "Read " . scalar(@names) . " species names from $infile" );
	}
	if ( $root_taxa ) {
		my @taxa = split(',', $root_taxa);
		push @names, @taxa;
		$log->info( scalar(@taxa) . " root taxa given in argument" );
	}	
	# expand root taxa if argument is provided
	if ( $expand_rank or $root_taxa ) {
		@names = $mts->expand_taxa( \@names, $expand_rank || "species" );
	}
	
	my @taxa_table = $mts->make_taxa_table( \@names, $opt->binomials_only );	
	$mts->write_taxa_file( $opt->outfile, @taxa_table );

	$log->info("DONE, results written to " . $opt->outfile);
	return 1;
}

1;
