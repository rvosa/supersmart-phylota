package Bio::SUPERSMART::App::smrtutils::Command::MarkerGraph;

use strict;
use warnings;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: writes dot language graph of marker table

=head1 NAME

MarkerGraph.pm - Writes dot language graph of marker table

=head1 SYNOPSYS

smrt-utils 

=head1 DESCRIPTION

Writes an undirected graph in dot (GraphViz) syntax that shows the connectedness
between taxa by way of the markers they share.

=cut

sub options {    
	my ($self, $opt, $args) = @_;
	my $infile_default = 'markers-backbone.tsv';
	return (
		['infile|i=s', "input marker table", { arg => 'file', default => $infile_default }],		
		['outfile|o=s', "output file in dot syntax", { arg => 'file' }],    	    
		['markers|m', "create an edge for each marker that connects a pair of taxa" ],
	    );	
}

sub validate {
	my ($self, $opt, $args) = @_;			

	$self->usage_error('need in file as argument') if not ($opt->infile);
	$self->usage_error('in file not found or empty') if not (-e $opt->infile and -s $opt->infile);
}

sub _openfh {
	my ( $log, $outfile ) = @_;
        my $outfh;
        if ( $outfile ) {
                open $outfh, '>', $outfile or die $!;
                $log->info("will write graph to outfile $outfile");
        }
        else {
                $outfh = \*STDOUT;
                $log->info("will write graph to STDOUT");
        }
	return $outfh;
}

sub _read_table {
	my ( $log, $infile ) = @_;
        $log->info("going to read marker table from $infile");
        my ( %marker, @header );
        open my $fh, '<', $infile or die $!;
        LINE: while(<$fh>) {
                chomp;
                if ( not @header ) {
                        @header = split /\t/, $_;
                        $marker{$header[$_]} = [] for 1 .. $#header;
                        next LINE;
                }
                my @record = split /\t/, $_;
                my $taxon = $record[0];
                for my $i ( 1 .. $#record ) {
                        if ( $record[$i] ) {
                                my $m = $header[$i];
                                push @{ $marker{$m} }, $taxon;
                                $log->debug("recorded marker $m for taxon $taxon");
                        }
                }
        }
	return %marker;
}

sub run {
	my ($self, $opt, $args) = @_;    
	my $log = $self->logger;
	my $infile = $opt->infile;
	my $outfile = $opt->outfile;
	my $markers = $opt->markers;	

	# read marker table
	my %marker = _read_table( $log, $infile );

	# open output handle
	my $outfh = _openfh( $log, $outfile );

	# write dot graph header
	print $outfh "graph markers {\n";

	# edges are individual markers
	if ( $markers ) {
        	for my $m ( keys %marker ) {
                	my @taxa = @{ $marker{$m} };
                	for my $i ( 0 .. $#taxa - 1 ) {
                        	for my $j ( $i + 1 .. $#taxa ) {
                                	print $outfh "\t", "\"$taxa[$i]\"", ' -- ', "\"$taxa[$j]\"", " [label=\"$m\"];\n";
                        	}
                	}
        	}
	}

	# edge weight is shared marker count
	else {
        	my %weight;
        	for my $m ( keys %marker ) {
                	my @taxa = @{ $marker{$m} };
                	for my $i ( 0 .. $#taxa - 1 ) {
                        	$weight{$taxa[$i]} = {} if not $weight{$taxa[$i]};
                        	for my $j ( $i + 1 .. $#taxa ) {
                                	$weight{$taxa[$i]}->{$taxa[$j]}++;
                        	}
                	}
        	}
        	for my $out_taxon ( keys %weight ) {
                	for my $in_taxon ( keys %{ $weight{$out_taxon} } ) {
                        	my $weight = $weight{$out_taxon}->{$in_taxon};
                        	print $outfh "\t", "\"$out_taxon\"", ' -- ', "\"$in_taxon\"", " [weight=$weight];\n";
                	}
        	} 
	}
	print $outfh "}\n";
}

1;
