package Bio::SUPERSMART::App::smrt::Command::BBdecompose;

use strict;
use warnings;

use List::MoreUtils 'uniq';
	
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector;
use Bio::Phylo::PhyLoTA::Service::TreeService;

use Bio::Phylo::PhyLoTA::Service::ParallelService;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
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
	my $outfile_default = "markers-clades.tsv";
	my $tree_default    = "consensus.nex";
	my $taxa_default    = "species.tsv";
	my $aln_default     = "aligned.txt";
	my $format_default  = "nexus";
	return (
		["backbone|b=s", "backbone tree as produced by 'smrt bbinfer'", { arg => "file", default => $tree_default}],
		["format|f=s", "file format of the backbone tree as produced by 'smrt consense'", { default => $format_default }],
		["classtree|c=s", "classification tree as produced by 'smrt classify', only needed with '-g' option", { arg => "file"}],
		["alnfile|a=s", "list of file locations of merged alignments as produced by 'smrt aln'", { arg => "file", default => $aln_default}],	
		["taxafile|t=s", "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'", { arg => "file", default => $taxa_default}],
		["add_outgroups|g", "attempt to automatically add outgroup from sister genus for each clade, if sufficient marker overlap between clades", { default=> 0} ],
		["outfile|o=s", "name of the output file (summary table with included accessions), defaults to $outfile_default", { default=> $outfile_default, arg => "file"}],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	#  If alignment or taxa file is absent or empty, abort  
	my @files = ( $opt->alnfile, $opt->backbone, $opt->taxafile );
	foreach my $file ( @files ){
		$self->usage_error("need alignment and taxa files and backbone tree file") if not $file;
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
	if ( $opt->add_outgroups ) {
		my $file = $opt->classtree;
		$self->usage_error("need classification tree if --add_outgroups option is set") if not $file;
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
	if ( $opt->format !~ /^(?:newick|nexus)$/i ) {
		$self->usage_error("only newick and nexus format are supported");
	}	
}

sub run{
	my ($self, $opt, $args) = @_;		
	
	# collect command-line arguments
	my $alnfile = $opt->alnfile;
	my $taxafile = $opt->taxafile;
	my $backbone = $opt->backbone;
	my $add_outgroup = $opt->add_outgroups;
	my $workdir = $self->workdir;
	my $outfile= $self->outfile;	

	# instantiate helper objects
	my $mt     = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
	my $mts    = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
	my $ts     = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $config = Bio::Phylo::PhyLoTA::Config->new;
	my $logger = $self->logger;
			
	# parse backbone tree
	$logger->info("going to read backbone tree $backbone");
	my $tree = parse_tree(
		'-format' => $opt->format,
		'-file'   => $backbone,
	);
	$ts->remap_to_ti($tree);
					
	# parse taxon mapping
	$logger->info("going to read taxa mapping $taxafile");
	my @records = $mt->parse_taxa_file($taxafile);
		
	# now read the list of alignments
	my @alignments;
	$logger->info("going to read list of alignments $alnfile");
	open my $fh, '<', $alnfile or die $!;
	while(<$fh>) {
		chomp;
		push @alignments, $_ if /\S/ and -e $_;
	}
	
	# decompose tree into clades and get the sets of species
	my @clades = map { my %h=('ingroup'=>$_); \%h; } $ts->extract_clades($tree, @records);;
	
	# get one outgroup species for each clade and append to species sets
	if ( $add_outgroup ){
		my $common = $opt->classtree;
		$logger->info("going to read classification tree $common");
		my $classtree = parse_tree(
			'-format'     => 'newick',
			'-file'       => $common,
			'-as_project' => 1,
		    );
		
		$ts->remap_to_ti($classtree);
		
		my $counter = 0;
		foreach my $clade ( @clades ){
			my $ingroup = $clade->{'ingroup'};
			my @og = $mts->get_outgroup_taxa( $classtree, $ingroup );
			
			# get the two species which occur in the most number of alignments
			my %aln_for_sp;
			for my $aln (@alignments) {
				my %fasta = $mt->parse_fasta_file($aln);
				for my $ogsp ( @og ) {
					my @sp = grep { /taxon\|$ogsp[^\d]/ } keys %fasta;
					@sp = map {$1 if $_=~m/taxon\|([0-9]+)/} @sp;
					$aln_for_sp{$_}++ for @sp;
				}
			}			
			my @sorted_sp = sort { $aln_for_sp{$a} <=> $aln_for_sp{$b} } keys %aln_for_sp;			
			my @outgroup = scalar @sorted_sp > 4 ? @sorted_sp[0..3] : @sorted_sp;

			$clade->{'outgroup'} = \@outgroup;
			$logger->info("Adding outgroup species " . join (', ', @outgroup) . " to clade #" . $counter);
			$counter++;
		}
	}

	# write suitable alignments to their respective clade folders
	#  and return a table with markers for all species
	my @table = pmap {
		my ($aln) = @_;
		my %fasta = $mt->parse_fasta_file($aln);		 	 
		$logger->info("checking whether alignment $aln can be included");

	        my $mdens = $config->CLADE_MIN_DENSITY;
		
		# return value: all included sequences per species
		my %ret;
		
		# assess for each set whether we have enough density
	      CLADE: for my $i ( 0 .. $#clades ) {

		      my %h = %{ $clades[$i] };
		      my %ingroup = map {$_=>1} @{$h{'ingroup'}};
		      my %outgroup = map {$_=>1} @{$h{'outgroup'}};
		      
		      my %seqs_ingroup = $mt->get_alignment_subset(\%fasta, {'taxon'=>[keys %ingroup]});		      
		      my %seqs_all = $mt->get_alignment_subset(\%fasta, {'taxon'=>[keys %outgroup, keys %ingroup]});		      
		      my $distinct = scalar ( keys %seqs_ingroup );


		      if ( $distinct  < 3 ) {
			      $logger->info("Not enough sequences in alignment $aln for clade # $i");
			      next CLADE;
		      }
		      # calculate distance of the subset of the alignment which includes only the (ingroup!) species in this clade		      		      
		      my $fastastr = '';
		      foreach my $k (keys %seqs_ingroup) {
			      $fastastr .= "> $k \n " . $seqs_ingroup{$k} . "\n";
		      }
		      my $dist  = $mt->calc_mean_distance($fastastr);
		      if ( $dist > $config->CLADE_MAX_DISTANCE ) {
			      $logger->info("$aln is too divergent for clade # $i, distance $dist > ".$config->CLADE_MAX_DISTANCE);
			      next CLADE;
		      }
		      		      
		      # the fraction of distinct, sequenced species is high enough,
		      # and the total number exceeds two (i.e. there is some topology to resolve)
		      if ( ($distinct/scalar(keys %ingroup)) >= $mdens ) {
			      $logger->info("$aln is informative and dense enough for clade $i");
			      # write alignment 
			      my ( $fh, $seed ) = _make_handle( $i, $aln, $workdir );
			      for my $defline ( keys %seqs_all ) {
				      print $fh '>', $defline, "\n", $fasta{$defline}, "\n";
				      my $species = $1 if $defline =~ /taxon\|([0-9]+)/;
				      $ret{$species} = [] if not $ret{$species};
				      push @{$ret{$species}}, $defline;					      				      
			      }
			      if ( $add_outgroup ) {
				      # write outgroup to file (skipped if already exists)
				      my @og = keys %outgroup;
				      _write_outgroup( $i, \@og, $workdir);
			      }
		      }
		      else {
			      my $dens = sprintf ( "%.2f", $distinct/scalar(keys %ingroup));
			      $logger->info("$aln is not informative or dense enough: density " 
					    . $dens . " < $mdens, number of distinct taxa in clade $i: $distinct");
		      }
	      }		
		return \%ret;
	} @alignments;
	
	# get all species that were present in an alignment and write marker table
	my @included_species = uniq map {keys(%$_)} @table;
	$mts->write_marker_summary( $outfile, \@table, \@included_species );

	$logger->info("DONE, results written into working directory $workdir");

	return 1;
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

# writes outgroup species into a file in the clade folder
sub _write_outgroup {
	my ( $i, $og, $workdir ) = @_;
	my $filename = "$workdir/clade$i/outgroup.txt";
	if  ( not -e $filename ) {
		open my $fh, '>', $filename or die $!;
		foreach my $sp ( @{$og} ) {
			print $fh  "$sp\n";
		}
		close $fh;
	}
}


1;
