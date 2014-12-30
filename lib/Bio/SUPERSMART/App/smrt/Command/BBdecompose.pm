package Bio::SUPERSMART::App::smrt::Command::BBdecompose;

use strict;
use warnings;

use List::MoreUtils 'uniq';
	
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Service::TreeService;

use Bio::Phylo::PhyLoTA::Service::ParallelService 'pthreads';

use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: decomposes a backbone tree into individual clades

=head1 NAME

BBdecompose.pm - decomposes a backbone tree into individual clades

=head1 SYNOPSYS

smrt bbdecompose [-h ] [-v ] [-w <dir>] -b <file> -c <file> -a <file> -t <file> [-g ] [-o <file>] 

=head1 DESCRIPTION

Given a rooted backbone phylogeny, a list of superclusters and a table of resolved taxa, 
decomposes the backbone into its constituent, most recent, monophyletic clades, expands 
these clades into all taxa from the provided table and assembles sets of alignments for 
each clade that can be used for further tree inference within that clade.

Traverses the backbone tree to find the nearest monophyletic clade that groups
the exemplar leaves. In the default case, the clade is the genus that subtends the
two exemplars provided they are monophyletic in the backbone tree. If they are not, we
traverse upwards to find the nearest monophyletic set of genera.

Each clade is then expanded into its constituent set of species on the basis of the
taxon mapping file taxa. For those sets of species, the list of alignment values
is evaluated, and for each alignment whose average divergence does not exceed 
CLADE_MAX_DISTANCE (given in the configuration file) but whose density of species 
in the sets does exceed CLADE_MIN_DENSITY
the sequences for the focal species are written to a new file, in a directory that
groups them with the other relevant alignments for that clade.

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "summary.tsv";
	return (
		["backbone|b=s", "backbone tree as produced by 'smrt bbinfer'", { arg => "file", mandatory => 1}],
		["classtree|c=s", "classification tree as produced by 'smrt classify'", { arg => "file", mandatory => 1}],
		["alnfile|a=s", "list of file locations of merged alignments as produced by 'smrt bbmerge'", { arg => "file", mandatory => 1}],	
		["taxafile|t=s", "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'", { arg => "file", mandatory => 1}],
		["add_outgroups|g", "automatically add outgroup for each clade", { default=> 0} ],
		["outfile|o=s", "name of the output file (summary table with included accessions), defaults to $outfile_default", { default=> $outfile_default, arg => "filename"}],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	#  If alignment or taxa file is absent or empty, abort  
	my @files = ( $opt->alnfile, $opt->backbone, $opt->taxafile, $opt->classtree );
	foreach my $file ( @files ){
		$self->usage_error("need alignment and taxa files and classification and backbone tree file") if not $file;
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
}

sub run{
	my ($self, $opt, $args) = @_;		
	
	# collect command-line arguments
	my $alnfile = $opt->alnfile;
	my $taxafile = $opt->taxafile;
	my $backbone = $opt->backbone;
	my $common = $opt->classtree;
	my $add_outgroup = $opt->add_outgroups;
	(my $workdir = $opt->workdir) =~ s/\/$//g;
	my $outfile= $self->outfile;	

	# instantiate helper objects
	my $mt = 'Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa';
	my $mts = 'Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector';
	my $ts = 'Bio::Phylo::PhyLoTA::Service::TreeService';
	my $config = Bio::Phylo::PhyLoTA::Config->new;
	my $logger = $self->logger;
			
	# parse backbone and NCBI common tree
	$logger->info("going to read backbone tree $backbone");
	my $tree = parse_tree(
		'-format'     => 'newick',
		'-file'       => $backbone,
		'-as_project' => 1,
	);
	$ts->remap_to_ti($tree);
		
	my $classtree = parse_tree(
		'-format'     => 'newick',
		'-file'       => $common,
		'-as_project' => 1,
	);
			
	$ts->remap_to_ti($classtree);
			
	# parse taxon mapping
	$logger->info("going to read taxa mapping $taxafile");
	my @records = $mt->parse_taxa_file($taxafile);
	
	# decompose tree into clades and get the sets of species
	my @set = $ts->extract_clades($tree, @records);
	
	# now read the list of alignments
	my @alignments;
	$logger->info("going to read list of alignments $alnfile");
	open my $fh, '<', $alnfile or die $!;
	while(<$fh>) {
		chomp;
		push @alignments, $_ if /\S/ and -e $_;
	}
	
	# get one outgroup species for each clade and append to species sets
	if ( $add_outgroup ){
		my $counter = 0;
		foreach (@set){
			my @og = $mts->get_outgroup_taxa( $classtree, $_ );
			push @{$_}, @og;
			$logger->info("Adding outgroup species " . join (', ', @og) . " to clade #" . ++$counter);
		}	
	}

	# write suitable alignments to their respective clade folders
	#  and return a table with markers for all species
	my @table = pmap {
	    my ($aln) = @_;
		$logger->info("checking whether alignment $aln can be included");
		my %fasta = $mt->parse_fasta_file($aln);
		my $dist  = $mt->calc_mean_distance(%fasta);
		my $nchar = $mt->get_nchar(%fasta);
		my $mdens = $config->CLADE_MIN_DENSITY;
		
		# return value: all included sequences per species
		my %ret;
		
		if ( $dist <= $config->CLADE_MAX_DISTANCE ) {
		
			# assess for each set whether we have enough density
			for my $i ( 0 .. $#set ) {
				my @species = @{ $set[$i] };
				my %seq;
				my $distinct = 0;
				for my $s ( @species ) {
					my @s = grep { /taxon\|$s[^\d]/ } keys %fasta;
					if ( @s ) {
						$seq{$s} = \@s;
						$distinct++;
					}
				}
				
				# the fraction of distinct, sequenced species is high enough,
				# and the total number exceeds two (i.e. there is some topology to resolve)
				if ( ($distinct/scalar(@species)) >= $mdens && $distinct > 2 ) {
					$logger->info("$aln is informative and dense enough for clade $i");
					my ( $fh, $seed ) = _make_handle( $i, $aln, $workdir );
					for my $s ( @species ) {
						if ( $seq{$s} ) {
							for my $j ( 0 .. $#{ $seq{$s} } ) {
								my $seq = $seq{$s}->[$j];
								print $fh '>', $seq, "\n", $fasta{$seq}, "\n";
								$ret{$s} = $seq;		
							}
						}
					}
				}
				else {
					my $dens = sprintf ( "%.2f", $distinct/scalar(@species));
					$logger->info("$aln is not informative or dense enough: density " 
					. $dens . " < $mdens, number of distinct taxa: $distinct");
				}
			}
		}
		else {
			$logger->info("$aln is too divergent: $dist > ".$config->CLADE_MAX_DISTANCE);
		}
		return \%ret;
	} @alignments;
	
	my @all_species = uniq map {@$_} @set;
	
	_write_marker_summary( $outfile, $mts, \@table, \@all_species );

	$logger->info("DONE, results written into working directory $workdir");

	return 1;
}


# writes a table containing all species as rows and all chosen markers 
#  as columns, reports the genbank accession of a specific marker for a species
sub _write_marker_summary {
	my ( $file, $mts, $tab, $specs ) = @_;
	my @table = @$tab;
	my @all_species = @$specs;
	
	# remove empty columns (for markers that are never included)
	@table = grep {keys %{$_}} @table;
	
	my @seed_gis = map { (values(%{$_}))[0] =~ /seed_gi\|([0-9]+)\|/ } @table;
	
	open my $outfh, '>', $file or die $!;
	
	# print table header
	print $outfh "taxon\t";
	print $outfh "marker" . $_ ."\t" for 1..$#seed_gis+1;
	print $outfh "\n";
	foreach my $species ( @all_species ) {
		my $name = $mts->find_node($species)->taxon_name;
		print $outfh $name . "\t";
		foreach my $column ( @table ) {
			my %h = %{$column};
			if (  my $seq = $h{$species} ) {
				# parse the GI from fasta descripion
				my ($gi) = $seq =~ /gi\|([0-9]+)\|/;
				my $seqobj = $mts->find_seq($gi);
				print $outfh $seqobj->acc . "\t";
			}
			else {
				print $outfh "\t";
			}
		}
		print $outfh "\n";
	}
	print $outfh "\n";
	
	# print information about markers at the bottom of the table
	foreach my $i ( 1..$#seed_gis+1 ) {
		my $seqobj = $mts->find_seq( $seed_gis[$i-1] );
		print $outfh "# marker$i cluster seed: " . $seqobj->acc . ", description: " . $seqobj->def . "\n"; 
	}
	close $outfh;
}

# makes file handle for alignment file to be written in clade directory
sub _make_handle {
	my ( $i, $aln, $workdir ) = @_;
	if ( not -d "$workdir/clade$i" ) {
		mkdir "$workdir/clade$i";
	}
	my ( $volume, $dir, $base ) = File::Spec->splitpath($aln);
	open my $fh, '>', "$workdir/clade$i/$base" or die $!;
	$base =~ s/\.fa$//;
	return $fh, $base;
}

1;
