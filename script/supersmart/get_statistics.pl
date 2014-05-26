#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use List::MoreUtils 'uniq';

use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::IO 'parse_tree';
use Bio::AlignIO;

=head1 NAME

get_statistics.pl - 

=head1 SYNOPSYS

 $ perl get_statistics.pl -t <taxa> -b <tree> -s <supermatrix> -w <dir> [--verbose]

=head1 DESCRIPTION

This script calculates basic statistics on the output of a SUPERSMART run.
Reported values include:

- number of distinct families in species list
- number of distinct genera in species list
- number of species
- number of tips in backbone tree
- % of gaps in backbone supermatrix
- number of supclades
- min and max alignment length for clade alignments
=cut

my $verbosity = WARN;
my ( $taxa, $supermatrix, $backbone, $workdir );
GetOptions(
	'taxa=s'     => \$taxa,
        'backbone=s' => \$backbone,
	'workdir=s'  => \$workdir,
        'supermatrix=s' => \$supermatrix,
	'verbose+'   => \$verbosity,
);


# instantiate helper objects
my $mt = 'Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa';
my $logger = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# get number of families, genera and species
my @records = $mt->parse_taxa_file($taxa);

my @levels = ("family", "genus", "species");

for my $level ( @levels ) {
        my @taxa = uniq ( map { $_->{$level} } @records );
        $logger->info("found ".scalar(@taxa)." distinct taxa of class $level");
}


# get number of terminals in backbone tree
my $tree = parse_tree(
	'-format'     => 'newick',
	'-file'       => $backbone,
	'-as_project' => 1,
);

my $ntax = $tree->get_ntax;
$logger->info("found $ntax tips in backbone tree");


# get statistics on backbone alignments. Note that here the supermatrix should be in FASTA format.
my $in  = Bio::AlignIO->new(-file   => $supermatrix ,
                            -format => 'fasta',); 
                            #-interleaved => 1,
                            #-idlength => 10);

#while ( my $aln = $in->next_aln() ) {
#        print "Alignment found!!\n";
#        print ref ($aln)."\n";
#}

my $align = $in->next_aln();

my $all_gaps = 0;
my $all_bp = 0;
foreach my $seq ($align->each_seq) {
        my $str = $seq->seq;
        my $count_mis = ($str =~ s/\?/\?/g);
        my $count_gap = ($str =~ s/\-/\-/g);
        my $count_A = ($str =~ s/A/A/g);
        my $count_C = ($str =~ s/C/C/g);
        my $count_T = ($str =~ s/T/T/g);
        my $count_G = ($str =~ s/G/G/g);
        my $count_tot= length($str) - $count_mis;
        my $perc_gap = $count_gap * 100 / $count_tot;
        $all_gaps += $count_gap;
        $all_bp += $count_tot;                
}

my $bb_perc = $all_gaps * 100 / $all_bp;
$logger->info("Backbone supermatrix : Total nuclotides : $all_bp total gaps : $all_gaps  % gaps : $bb_perc");

# count number of subclades
opendir (DIR, $workdir) or die ("Could not open working directory $workdir");
my @cladedirs = grep m/^clade/,  readdir DIR;
close DIR;

my $clade_count = scalar @cladedirs;
$logger->info("Number of subclades : $clade_count");

# get min and max base pairs of all alignments for all clades
my $min_seqlength = 1000000;
my $max_seqlength = -1;
foreach my $dir (@cladedirs){
        opendir (CLDIR, $workdir."/".$dir) or die ("Could not open clade directory");
        my @aln_files = grep m/\.fa/, readdir CLDIR;
        close CLDIR;
        foreach my $file (@aln_files) {
                my $in  = Bio::AlignIO->new(-file   => $workdir."/".$dir."/".$file ,
                                             -format => 'fasta',); 
                my $aln = $in->next_aln();
                foreach my $seq($aln->each_seq){
                        my $length = length($seq->seq);
                        if ($length <= $min_seqlength) {
                                $min_seqlength = $length;
                        }
                        if ($length >= $max_seqlength) {
                                $max_seqlength = $length;
                        }
                }
                
        }
}

$logger->info("Min length of all clade alignments : $min_seqlength");
$logger->info("Max length of all clade alignments : $max_seqlength");
