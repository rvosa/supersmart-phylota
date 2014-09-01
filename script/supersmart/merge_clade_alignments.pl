#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::PhyLoTA::Service;
use Bio::Phylo::IO 'parse_matrix';
use Bio::Phylo::Util::Logger ':levels';

=head1 NAME

merge_clade_alignments.pl - creates a merged data file for species-level inference

=head1 SYNOPSYS

 $ perl merge_clade_alignments.pl -w <dir> [--verbose]

=head1 DESCRIPTION

Given a working directory, traverses it looking for subdirectories of the pattern
C<clade\d+>. Perusing each of these, it merges the *.fa (FASTA) files in it and
produces a single output file that can be analysed by the infer_clade.pl wrapper.

=cut

my $ns = 'http://www.supersmart-project.org/terms#';

# process command line arguments
my $verbosity = WARN;
my ( $workdir );
GetOptions(
	'verbose+' => \$verbosity,
	'workdir=s'=> \$workdir,
);

# instantiate helper objects
my $factory = Bio::Phylo::Factory->new;
my $service = Bio::Phylo::PhyLoTA::Service->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# start iterating
$log->info("going to look for clade data in $workdir");
opendir my $odh, $workdir or die $!;
while( my $dir = readdir $odh ) {
	if ( $dir =~ /^clade\d+$/ and -d "${workdir}/${dir}" ) {
	
		# start processing the directory
		$log->info("going to merge alignments in $workdir/$dir");
		my $project = $factory->create_project( '-namespaces' => { 'smrt' => $ns } );
		my $taxa = $factory->create_taxa;
		$project->insert($taxa);
		my %t; # to keep track of taxa
		opendir my $idh, "${workdir}/${dir}" or die $!;
		while( my $file = readdir $idh ) {
		
			# found an alignment
			if ( $file =~ /(.+?)\.fa$/ ) {
				my $id = $1;
				
				# parse the file
				$log->info("adding alignment $file");
				my $matrix = parse_matrix(
					'-type'       => 'dna',
					'-format'     => 'fasta',
					'-file'       => "${workdir}/${dir}/${file}",
					'-as_project' => 1,
				);
				$matrix->set_name($id);
												
				# update sequence labels, link to taxon objects
				my ( %gaps, $ntax );
				$matrix->visit(sub{
					my $row = shift;
					my $name = $row->get_name;
					my %fields = split /\|/, $name;
					my $binomial = $service->find_node($fields{'taxon'})->taxon_name;
					$binomial =~ s/ /_/g;
					
					# create new taxon object if none seen yet
					if ( not $t{$binomial} ) {
						$t{$binomial} = $factory->create_taxon( '-name' => $binomial );
						$t{$binomial}->set_meta_object( 'smrt:tid' => $fields{'taxon'} );
						$taxa->insert($t{$binomial});
					}
					$row->set_name( $binomial );
					$row->set_taxon( $t{$binomial} );
					$row->set_meta_object( 'smrt:gi' => $fields{'gi'} );
					
					# keep track of which columns might be all gaps
					my @char = $row->get_char;
					for my $i ( 0 .. $#char ) {
						my $c = $char[$i];
						if ( $c eq '-' or $c eq '?' or $c eq 'N' ) {
							$gaps{$i}++;							
						}
					}
					$ntax++;
				});
				$matrix->set_taxa($taxa);
				
				# delete all-gaps columns. this can happen because we are now 
				# operating on taxon subsets of the clustered alignments
				my $max = $matrix->get_nchar - 1;
				my @indices = grep { not defined $gaps{$_} or $gaps{$_} < $ntax } 0 .. $max;
				$matrix->visit(sub{
					my $row = shift;
					my @char = $row->get_char;
					my @ungapped = @char[@indices];
					$row->set_char(@ungapped);
				});
				
				# add to project
				$project->insert($matrix);
			}
		}
		
		# write the merged nexml
		$log->info("going to write file ${workdir}/${dir}/${dir}.xml");
		open my $outfh, '>', "${workdir}/${dir}/${dir}.xml" or die $!;
		print $outfh $project->to_xml( '-compact' => 1 );
	}
}
