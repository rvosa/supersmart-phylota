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
        ["backbone|b=s", "backbone tree as produced by e.g. 'smart bbinfer' and 'smrt consense', defaults to $tree_default", { arg => "file", default => $tree_default}],
        ["format|f=s", "file format of the backbone tree as produced by 'smrt consense', defaults to $format_default", { default => $format_default, arg => "format" }],
        ["classtree|c=s", "classification tree as produced by 'smrt classify', only needed with '-g' option", { arg => "file"}],
        ["alnfile|a=s", "list of file locations of merged alignments as produced by 'smrt aln'", { arg => "file", default => $aln_default}],    
        ["taxafile|t=s", "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'", { arg => "file", default => $taxa_default}],
        ["add_outgroups|g", "attempt to automatically add outgroup from sister genus for each clade, if sufficient marker overlap between clades", { default=> 0 } ],
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
    my $alnfile      = $opt->alnfile;
    my $taxafile     = $opt->taxafile;
    my $backbone     = $opt->backbone;
    my $add_outgroup = $opt->add_outgroups;
    my $common       = $opt->classtree;
    my $workdir      = $self->workdir;
    my $outfile      = $self->outfile;  

    # instantiate helper objects
    my $mt     = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
    my $mts    = Bio::Phylo::PhyLoTA::Service::MarkersAndTaxaSelector->new;
    my $ts     = Bio::Phylo::PhyLoTA::Service::TreeService->new;
    my $config = Bio::Phylo::PhyLoTA::Config->new;
    my $logger = $self->logger;
            
    # parse backbone tree
    $logger->info("Going to read backbone tree $backbone");
    my $tree = parse_tree(
        '-format' => $opt->format,
        '-file'   => $backbone,
    );
    $ts->remap_to_ti($tree);
                    
    # parse taxon mapping
    $logger->info("Going to read taxa mapping $taxafile");
    my @taxa = $mt->parse_taxa_file($taxafile);
        
    # now read the list of alignments
    my @alignments;
    $logger->info("Going to read list of alignments $alnfile");
    open my $fh, '<', $alnfile or die $!;
    while(<$fh>) {
        chomp;
        push @alignments, $_ if /\S/ and -e $_;
    }
    
    # decompose tree into clades and get the sets of species
    my @clades = map { { 'ingroup' => $_ } } $ts->extract_clades($tree, @taxa);
    
    # get the exemplars
    for my $c ( @clades ) {
        my @ex = grep { defined $_ } map { $tree->get_by_name($_) } @{ $c->{'ingroup'} };
        $c->{'exemplars'} = [ keys %{{ map { $_->get_name => 1 } @ex }} ];
    }
 
    # get one outgroup species for each clade and append to species sets
    if ( $add_outgroup ) {
        
        # read the classification tree, source of candidate outgroups
        $logger->info("Going to read classification tree $common");
        my $classtree = parse_tree(
            '-format'     => 'newick',
            '-file'       => $common,
            '-as_project' => 1,
        );      
        $ts->remap_to_ti( $classtree, @taxa );
        
        # iterate over hashes, with key 'ingroup', value is an
        # array ref of taxon IDs
        my $counter = 0;
        for my $clade ( @clades ){
            my $ingroup = $clade->{'ingroup'};
            my @og = $mts->get_outgroup_taxa( $classtree, $ingroup );
            
            # get the two species which occur in the most number of alignments
            my %aln_for_sp;
            for my $aln ( @alignments ) {
                my %fasta = $mt->parse_fasta_file($aln);
                for my $ogsp ( @og ) {
                    my @sp = grep { /taxon\|$ogsp[^\d]/ } keys %fasta;
                    @sp = map { $1 if $_ =~ m/taxon\|([0-9]+)/ } @sp;
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
	
	# collect alignments for clades and write them into respective clade directories
	my @table = pmap {
		
		my $i = $_;

		# for each clade, collect a set of alignments   
	    my $clade = $clades[$i];
	    my %ingroup  = map { $_ => 1 } @{ $clade->{'ingroup'} };
		my %outgroup = map { $_ => 1 } @{ $clade->{'outgroup'} };  
		my %exemplars = map { $_ => 1 } @{ $clade->{'exemplars'} };  

		# write outgroup to file (skipped if already exists)
		_write_outgroup($i,[keys %outgroup],$workdir) if $add_outgroup;

	    my $mindens = $config->CLADE_MIN_DENSITY;
	    my $maxdist = $config->CLADE_MAX_DISTANCE;	    
	    
	    # loop over alignments and assess whether it is suitable for the clade
	    my @clade_alignments;
	    ALN: for my $aln (@alignments) {
		    
		    # make subset: take only the sequences that are in the clade (or outgroup, if given)
		    my %fasta = $mt->parse_fasta_file($aln);             
		    $logger->info("Checking whether alignment $aln can be included");
		   
		    my %seqs_ingroup = $mt->get_alignment_subset(\%fasta, {'taxon'=>[keys %ingroup]});            
		    my %seqs_all = $mt->get_alignment_subset(\%fasta, {'taxon'=>[keys %outgroup, keys %ingroup]});            
		    
		    # check if density is high enough
		    my $distinct = scalar keys %seqs_ingroup;		    
		    if ( ($distinct/scalar keys %ingroup) < $mindens ) {
			    my $dens = sprintf "%.2f", $distinct / scalar keys %ingroup;
			    $logger->info("$aln is not  dense enough (density " . $dens . " < $mindens) for clade # $i");
			    next ALN;
		    }
		    
		    # check if distance is not too high
		    my $dist = $mt->calc_mean_distance($mt->to_fasta_string(%seqs_ingroup));
		    if ( $dist > $maxdist ) {
			    $logger->info("$aln is too divergent (distance $dist > $maxdist) for clade # $i");
			    next ALN;
		    }
		    
		    # add alignment to set of clade alignments
		    push @clade_alignments, \%seqs_all;		   
	    }
		
	    # for the given set of taxa and alignments, get all subsets of taxa that share at least one marker
		my @all_taxa = map {@$_} map {values(%$_)} @clades;
	    my %adj = $mts->generate_marker_adjacency_matrix(\@clade_alignments, [keys %outgroup, keys %ingroup]);
		my @subsets = @{$mts->get_connected_taxa_subsets(\%adj)};
	
		# sort subsets by size, decreasing
		@subsets = sort { scalar(@$b) <=> scalar(@$a)  } @subsets;
		
		# select the largest subset with exemplars present
		my @set;
		for my $s (@subsets) {
			if ( scalar grep( $exemplars{$_}, @$s)  == 2) {
				@set = @$s;
				last;
			}
		}
		if ( ! scalar (@set) ) {
			$logger->warn("Could not find sufficient data for exemplar species. Skipping clade with taxa " . join(',', keys(%ingroup)));
		}
		else {
			my @table_for_aln = _write_clade_alignments( $i, \@clade_alignments, \@set, $self->workdir );
			return @table_for_aln;
		}
		
	}  ( 0..$#clades);
    
    # get all species that were present in an alignment and write marker table
    my @included_species = uniq map {keys(%$_)} @table;
	@included_species = sort(@included_species);
    $mts->write_marker_table( $outfile, \@table, \@included_species );

    $logger->info("DONE, results written into working directory $workdir");

    return 1;
}

# writes alignments for a given set of taxa
sub _write_clade_alignments {
	my ( $clade, $alns, $taxa, $workdir ) = @_;

	# will return a markers table; rows are taxa, columns are different markers
	my @table;
	
	if ( not -d "$workdir/clade$clade" ) {
        mkdir "$workdir/clade$clade";
    }

	my @alignments = @$alns;
	for my $i (0..$#alignments) {
		my %aln = %{$alignments[$i]};
		my $idx = $i+1;
		open my $fh, '>', "$workdir/clade$clade/clade$clade-aln$idx.fa" or die $!;
		# column for markers table
		my %column;
		for my $defline( keys %aln ) {
			my $species = $1 if $defline =~ /taxon\|([0-9]+)/;
			
            # only write sequence if taxon is in the given list of taxa
			if ( grep($species, @$taxa) ) {
				print $fh '>', $defline, "\n", $aln{$defline}, "\n";
			}
			$column{$species} = [] if not $column{$species};
			push @{ $column{$species} }, $defline;
		}
		push @table, \%column;
		close $fh;
	}
	return @table;		
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
