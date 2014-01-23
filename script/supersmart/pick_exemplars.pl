#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my ( $taxatable, $fastafile );
GetOptions(
    'verbose+'    => \$verbosity,
    'taxatable=s' => \$taxatable,
    'fastafile=s' => \$fastafile,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
    '-class' => 'main',
    '-level' => $verbosity,
);

# read taxa file
my %genus_for;
{
    $log->info("going to read taxa mapping from $taxatable");
    open my $fh, '<', $taxatable or die $!;
    my @header;
    LINE: while(<$fh>) {
        chomp;
        my @line = split /\t/, $_;
        if ( not @header ) {
            @header = @line;
            next LINE;
        }
        my ( $species, $genus );
        for my $i ( 0 .. $#header ) {
            if ( $header[$i] eq 'species' ) {
                $species = $line[$i];
            }
            if ( $header[$i] eq 'genus' ) {
                $genus = $line[$i];
            }
        }
        $genus_for{$species} = $genus;
    }
    $log->info("done reading taxa mapping");
}

# read fasta
my %fasta;
{
    $log->info("going to read sequences from $fastafile");
    open my $fh, '<', $fastafile or die $!;
    %fasta = simple_fasta( do { local $/; <$fh> } );
    close $fh;
    $log->info("done reading sequences");
}

# cluster def lines by genus
my %seqs_for_genus;
{
    $log->info("going to cluster definition lines by genus");
    for my $defline ( keys %fasta ) {
        if ( $defline =~ /taxon\|(\d+)/ ) {
            my $species = $1;
            if ( my $genus = $genus_for{$species} ) {
                $seqs_for_genus{$genus} = [] if not $seqs_for_genus{$genus};
                push @{ $seqs_for_genus{$genus} }, $defline;
            }
        }
    }
}

# iterate over genera, write exemplar
$log->info("going to pick exemplars");
for my $genus ( keys %seqs_for_genus ) {
    my $length = 0;
    my $longest;
    for my $defline ( @{ $seqs_for_genus{$genus} } ) {
        my $seq = $fasta{$defline};
        $seq =~ s/-//g;
	$seq =~ s/\?//g;
        if ( length($seq) > $length ) {
            $longest = $defline;
            $length  = length($seq);
        }
    }
    print '>', $longest, "\n", $fasta{$longest}, "\n";
}

# reads a FASTA string, returns a hash keyed on definition line (sans '>'
# prefix), value is concatenated seq
sub simple_fasta {
    my $string = shift;
    my @lines = split /\n/, $string;
    my %fasta;
    my $current;
    for my $line ( @lines ) {
        chomp $line;
        if ( $line =~ /^>(.+)/ ) {
            $current = $1;
            if ( exists $fasta{$current} ) {
                $log->warn("already seen definition line $current");
                $fasta{$current} = '';
            }            
        }
        else {
            $fasta{$current} .= $line;
        }
    }
    return %fasta;
}
