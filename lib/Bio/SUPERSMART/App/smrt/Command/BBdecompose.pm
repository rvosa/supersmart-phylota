package Bio::SUPERSMART::App::smrt::Command::BBdecompose;

use strict;
use warnings;

use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use Bio::Phylo::PhyLoTA::Service::ParallelService 'pthreads'; # can be either 'pthreads' or 'mpi';


use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);



sub options {
	my ($self, $opt, $args) = @_;		
	return (
		["backbone|b=s", "backbone tree as produced by 'smrt bbinfer'", { arg => "file", mandatory => 1}],
		["commontree|c=s", "classification tree as produced by 'smrt classify'", { arg => "file", mandatory => 1}],
		["alnfile|a=s", "list of file locations of merged alignments  as produced by 'smrt orthologize'", { arg => "file", mandatory => 1}],	
		["taxafile|t=s", "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'", { arg => "file", mandatory => 1}],
		[ "add_outgroups|o", "automatically add outgroup for each clade", { default=> 1} ],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	#  If alignment or taxa file is absent or empty, abort  
	my @files = ( $opt->alnfile, $opt->backbone, $opt->taxafile, $opt->commontree );
	foreach my $file ( @files ){
		print $file . "\n";
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
	my $common = $opt->commontree;
	my $add_outgroup = $opt->add_outgroups;
	(my $workdir = $opt->workdir) =~ s/\/$//g;
	#my $outfile= $self->outfile;	

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
	my $commontree = parse_tree(
		'-format'     => 'newick',
		'-file'       => $common,
		'-as_project' => 1,
	);
			
			
	# parse taxon mapping
	$logger->info("going to read taxa mapping $taxafile");
	my @records = $mt->parse_taxa_file($taxafile);
	
	# decompose tree into clades and get the sets of species
	my @set = $ts->extract_clades($tree, @records);
	
	# now read the list of alignments
	my @alignments;
	sequential {
		$logger->info("going to read list of alignments $alnfile");
		open my $fh, '<', $alnfile or die $!;
		while(<$fh>) {
			chomp;
			push @alignments, $_ if /\S/ and -e $_;
		}
	};
	
	# get one outgroup species for each clade and append to species sets
	if ( $add_outgroup ){
		my $counter = 0;
		foreach (@set){
			my @og = $mts->get_outgroup_taxa( $commontree, $_ );
			push @{$_}, @og;
			$logger->info("Adding outgroup species " . join (' ,', @og) . " to clade #" . ++$counter);
		}	
	}
	
	# write suitable alignments to their respective clade folders
	pmap {
	    my $aln = $_;
		$logger->info("assessing alignment file $aln");
		my %fasta = $mt->parse_fasta_file($aln);
		my $dist  = $mt->calc_mean_distance(%fasta);
		my $nchar = $mt->get_nchar(%fasta);
		my $mdens = $config->CLADE_MIN_DENSITY;
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
							}
						}
					}
				}
			}
		}
		else {
			$logger->info("$aln is too divergent: $dist > ".$config->CLADE_MAX_DISTANCE);
		}
		$logger->info("done assessing alignment file $aln");
	} @alignments;
	
	$logger->info("DONE, results written to $workdir");

	return 1;
}

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
