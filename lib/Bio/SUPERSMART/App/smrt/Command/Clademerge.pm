package Bio::SUPERSMART::App::smrt::Command::Clademerge;

use strict;
use warnings;

use Bio::Phylo::Factory;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::IO 'parse_matrix';
use Bio::Phylo::Util::CONSTANT ':objecttypes';

use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: merges sets of alignments into input files for clade inference

=head1 NAME

Clademerge.pm - For each decomposed clade, merges the set of alignments assembled for this clade into an input file for tree inference.

=head1 SYNOPSYS

smrt clademerge [-h ] [-v ] [-w <dir>] [-o <format>]

=head1 DESCRIPTION

Given a working directory, traverses it looking for subdirectories of the pattern
clade*. Perusing each of these, it merges the *.fa (FASTA) files in it and
produces a single output file that can be analysed by the subcommand cladeinfer.

=cut


sub options {
	my ($self, $opt, $args) = @_;		
	return (
		[ "outformat|o=s", "output format for merged clade files (phylip or nexml), defaults to 'nexml'", { arg=>"format", default=> 'nexml'} ],
		[ "outfile|f=s", "location of output directory", {arg=>"location", default => "clademerge_out.txt"} ],
	);	
}

sub validate {};

sub run {
	my ($self, $opt, $args) = @_;		
	
	my $workdir = $self->workdir;
	my $outformat = $opt->outformat;
	my $outfile= $self->outfile;	
		
	# instantiate helper objects
	my $factory = Bio::Phylo::Factory->new;
	my $service = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $log = $self->logger;
	my $ns = 'http://www.supersmart-project.org/terms#';

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
						my $binomial = $fields{'taxon'}; ##$service->find_node($fields{'taxon'})->taxon_name;
						#$binomial =~ s/ /_/g;
						
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
			
			if (lc $outformat eq 'nexml'){
				# write the merged nexml
				$log->info("going to write file ${workdir}/${dir}/${dir}.xml");
				open my $outfh, '>', "${workdir}/${dir}/${dir}.xml" or die $!;
				print $outfh $project->to_xml( '-compact' => 1 );
			}
			elsif (lc $outformat eq 'phylip'){
				my @matrices = @{ $project->get_items(_MATRIX_) };
				my ($taxa) = @{$project->get_items(_TAXA_)} ;
				$log->info("going to write file ${workdir}/${dir}/${dir}.phy");
				$service->make_phylip_from_matrix($taxa, "${workdir}/${dir}/${dir}.phy", @matrices);
			}
		}
	}
	open my $fh, '>', $outfile;
	print $fh "Clademerge done\n";
	close $fh;
	$log->info("DONE");
	return 1;
}

1;