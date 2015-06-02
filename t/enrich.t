#!/usr/bin/perl
use strict;
use warnings;
use Test::More 'no_plan';
use FindBin '$Bin';
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse_matrix';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;

my $outfile = shift;

my $mts     = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
my $fac     = Bio::Phylo::Factory->new;
my $ns      = 'http://www.supersmart-project.org/terms#';
my $file    = "${Bin}/testdata/186701483-16739-4-subtree.fa";
my $project = $fac->create_project( '-namespaces' => { 'smrt' => $ns } );
my $taxa    = $fac->create_taxa;
my $matrix  = parse_matrix(
	'-type'       => 'dna',
        '-format'     => 'fasta',
        '-file'       => $file,
        '-as_project' => 1,
);

$project->insert($taxa);
$project->insert($matrix);

my %t;
$matrix->visit(sub{
	my $row = shift;
        my $name = $row->get_name;
        my %fields = split /\|/, $name;
        $fields{$_} =~ s/^(\d+).*$/$1/ for keys %fields;
        my $binomial = $fields{'taxon'};

        # create new taxon object if none seen yet
        if ( not $t{$binomial} ) {
        	$t{$binomial} = $fac->create_taxon( '-name' => $binomial );
                $t{$binomial}->set_meta_object( 'smrt:tid' => $fields{'taxon'} );
                $taxa->insert($t{$binomial});
        }
        $row->set_name( $binomial );
        $row->set_taxon( $t{$binomial} );
        $row->set_meta_object( 'smrt:gi'      => $fields{'gi'} );
        $row->set_meta_object( 'smrt:mrca'    => $fields{'mrca'} );
        $row->set_meta_object( 'smrt:seed_gi' => $fields{'seed_gi'} );
});
$matrix->set_taxa($taxa);

$mts->enrich_matrix($matrix);

my $xml = $project->to_xml( '-compact' => 1 );

ok( $xml );

if ( $outfile ) {
	open my $fh, '>', $outfile or die $!;
	print $fh $xml;
	close $fh;
}
