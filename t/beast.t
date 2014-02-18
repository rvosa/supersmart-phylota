#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';
use Bio::Phylo::IO 'parse_matrix';
use Bio::Phylo::Project;

# test instantiation
BEGIN { use_ok('Bio::Tools::Run::Phylo::StarBEAST'); }
my $beast = new_ok('Bio::Tools::Run::Phylo::StarBEAST');

# make list of alignments
my $p = Bio::Phylo::Project->new;
opendir my $dh, "${Bin}/beast" or die $!;
while( my $entry = readdir $dh ) {
	next unless $entry =~ /\.fa$/;
	my $matrix = parse_matrix(
		'-format' => 'fasta',
		'-file'   => "${Bin}/beast/$entry",
		'-type'   => 'dna',
	);
	$entry =~ s/\.fa$//;
	$matrix->set_name($entry);
	$p->insert($matrix);
}
ok($beast->_alignment($p), "read 7 alignments");

my $twig = $beast->_make_beast_xml;
$twig->print;