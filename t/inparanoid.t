#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::SequenceGetter;

# make logger more verbose during run_blast_search
my $log = Bio::Phylo::Util::Logger->new(
	'-level'  => DEBUG,
	'-method' => 'Bio::Phylo::PhyLoTA::Service::SequenceGetter::run_blast_search',
);

# instantiate sequence service object
my $sg = Bio::Phylo::PhyLoTA::Service::SequenceGetter->new;

# AAEL004420-PB__spec_id_1
my $seq = 'MDSRKNKRSSNSKQKNASKQDSKSITNKPLPDGAGSQPRILLKPKENAEVLRKPEPVITIATKAPP' .
'AKPDPSHHAAKDTAVLKQRPVAQEIPQSITIVNSMLKSVNLVNGNNALNAHAFDYLHESNSDFFVIGAIGMQGSGK' .
'STVLNLLAADPTEEAIRQAAFHVGGGVFPITSLFGRSETVEFEDAEIRMHITKDRIILLDSASVLSNRGNKDFVLS' .
'ELDDIRRIMLLLSVCHVLLVLQEDYFNINFIRLLRCAEMMIQRDQKDTQHLNPRIIFVKNKCNRNNFTIQDKTFHE' .
'KIYKQILKDTKLKIYSNEDDKEKINAIYLPKLFYDNLLFTDDDSSFSCILKLRQHVFMTPTHDAVEGCESLTEKAW' .
'SQIVHHVLESHHNNYFLRKYENLKEKYNLHNHVNVVENAAKEKSYLNFIDT';

# run blast
my @hits = $sg->run_blast_search( '-seq' => $seq );

ok( @hits );