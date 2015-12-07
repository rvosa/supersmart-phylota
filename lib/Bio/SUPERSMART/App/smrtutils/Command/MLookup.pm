package Bio::SUPERSMART::App::smrtutils::Command::MLookup;
use strict;
use warnings;
use Template;
use File::Spec;
use List::Util 'sum';
use List::MoreUtils qw(pairwise uniq all);
use Bio::Phylo::Treedrawer;
use Bio::Phylo::IO qw(parse_tree);
use Bio::SUPERSMART::Service::DecorationService;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);


# ABSTRACT: marker lookup

=head1 NAME

MLookup - Looks up marker names for accession numbers

=head1 SYNOPSYS

=head1 DESCRIPTION

Looks up marker names (e.g. gene names or other short descriptors)
for accession numbers as listed in the backbone markers table, the
clade markers table, or on the command line.

=cut

sub options {
	my ($self, $opt, $args) = @_;
	my $bbmarkers_default = 'markers-backbone.tsv';
	my $clade_default     = 'markers-clades.tsv';
	my $outfile_default   = 'markers.tsv';
	return (
		['backbone|b=s', "backbone markers table. Default: $bbmarkers_default", { default => $bbmarkers_default, arg => 'file' } ],
		['clades|c=s', "clade markers table. Default: $clade_default", { default => $clade_default, arg => 'file' } ],
		['outfile|o=s', "output file. Default: $outfile_default", { default => $outfile_default, arg => 'file' } ],
		['accession|a=s', "accession numbers, comma-separated", { arg => 'acc1,acc2,...' } ],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		
	
	# check if exists
	if ( not -e $opt->backbone ) {
		$self->logger->info("Backbone marker table not found, will not process ".$opt->markers); 
	}
	if ( not -e $opt->clades ) {
		$self->logger->info("Clade marker table not found, will not process ".$opt->clades);
	}
}

sub run {
	my ($self, $opt, $args) = @_;
	my $outfile = $opt->outfile;
	my $logger  = $self->logger;	

	# collect all distinct accession numbers	
	my ( %acc, @acc );		
	if ( -e $opt->backbone and -s $opt->backbone ) {
		$logger->info("Going to lookup backbone markers from file ".$opt->backbone);
		push @acc, $self->_get_acc($opt->backbone);
	}
        if ( -e $opt->clades and -s $opt->clades ) {
                $logger->info("Going to lookup clade markers from file ".$opt->clades);
		push @acc, $self->_get_acc($opt->clades);
        }
	if ( $opt->accession ) {
		$logger->info("Going to lookup markers from command line: ".$opt->accession);
		push @acc, split /,/, $opt->accession;
	}
	%acc = map { $_ => 0 } @acc;
	$logger->info("Going to lookup ".scalar(keys(%acc))." accession numbers");

	# do the lookup
	my $sg = Bio::SUPERSMART::Service::SequenceGetter->new;
	for my $acc ( keys %acc ) {
		my @result = $sg->get_markers_for_accession($acc);
		$acc{$acc} = join ', ', @result;
		$logger->info("$acc => $acc{$acc}");
	}

	# write to file
	open my $fh, '>', $outfile or die $!;
	for my $acc ( keys %acc ) {
		print $fh $acc, "\t", $acc{$acc}, "\n";
	}
	close $fh;
	
	$logger->info("DONE. Outfile written to $outfile");
	return 1;
}

sub _get_acc {
	my ( $self, $file ) = @_;
        my $mt = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
	my @records = $mt->parse_taxa_file($file);
	my @result;
	for my $r ( @records ) {
		for my $key ( @{ $r->{'keys'} } ) {
			next if $key eq 'taxon';
			if ( $r->{$key} ) {
				push @result, $r->{$key};
			}
		}
	}
	return @result;
}

1;
