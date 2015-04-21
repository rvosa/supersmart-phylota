package Bio::Tools::Run::Phylo::StarBEAST;
use strict;
use version;
use XML::Twig;
use File::Temp 'tempfile';
use Bio::AlignIO;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Forest::Tree;
use Bio::Align::DNAStatistics;
use Bio::Tree::DistanceFactory;
use Bio::Tools::Run::Phylo::PhyloBase;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use base qw(Bio::Tools::Run::Phylo::PhyloBase);

our $PROGRAM_NAME = 'beast';
our @beast_PARAMS = (
	'mc3_chains',       # number of chains
	'mc3_delta',        # temperature increment parameter
	'mc3_temperatures', # a comma-separated list of the hot chain temperatures
	'mc3_swap',         # frequency at which chains temperatures will be swapped
	'seed',             # Specify a random number generator seed
	'threshold',        # Full evaluation test threshold (default 1E-6)
	'beagle_order',     # set order of resource use
	'beagle_scaling',   # specify scaling scheme to use
	'beagle_rescale',   # frequency of rescaling (dynamic scaling only)	
	'beagle_instances', # divide site patterns amongst instances	
);
our @beast_SWITCHES = (
	'verbose',          # Give verbose XML parsing messages
	'warnings',         # Show warning messages about BEAST XML file
	'strict',           # Fail on non-conforming BEAST XML file
	'beagle_CPU',       # use CPU instance		
	'beagle_GPU',       # use GPU instance if available
	'beagle_SSE',       # use SSE extensions if available
	'beagle_cuda',      # use CUDA parallization if available
	'beagle_opencl',    # use OpenCL parallization if available
	'beagle_single',    # use single precision if available
	'beagle_double',    # use double precision if available	
	'overwrite',        # Allow overwriting of log files
);
my $beast_ns  = 'http://beast.bio.ed.ac.uk/BEAST_XML_Reference#';
my $beast_pre = 'beast';
my $fac = Bio::Phylo::Factory->new;
my $log = Bio::Phylo::Util::Logger->new;

=head1 NAME

Bio::Tools::Run::Phylo::StarBEAST - module for interface with the BEAST software 

=head1 SYNOPSYS

 use Bio::Tools::Run::Phylo::StarBEAST;
 use FindBin '$Bin';
 
 my $beast = Bio::Tools::Run::Phylo::StarBEAST->new;
 $beast->chain_length(100);
 $beast->sample_freq(100);
 $beast->run( "$Bin/beast.xml" );


=head1 DESCRIPTION

This module provides an interface to the  *BEAST software
(http://beast.bio.ed.ac.uk) for inference of phylogenetic trees 
using Markov-Chain Monte-Carlo methods.   

=head1 METHODS

In addition to the methods described below, all named arguments described for the
constructor are available as methods as well. For example, the argument C<-verbose>
becomes:

 $beast->verbose(1);

=over

=item new

The constructor takes the following optional named arguments that require a numerical 
value:

	'-mc3_chains',       # number of chains
	'-mc3_delta',        # temperature increment parameter
	'-mc3_temperatures', # a comma-separated list of the hot chain temperatures
	'-mc3_swap',         # frequency at which chains temperatures will be swapped
	'-seed',             # Specify a random number generator seed
	'-threshold',        # Full evaluation test threshold (default 1E-6)
	'-beagle_order',     # set order of resource use
	'-beagle_scaling',   # specify scaling scheme to use
	'-beagle_rescale',   # frequency of rescaling (dynamic scaling only)	
	'-beagle_instances', # divide site patterns amongst instances	

In addition, the following named arguments that require a boolean are available:

	'-verbose',          # Give verbose XML parsing messages
	'-warnings',         # Show warning messages about BEAST XML file
	'-strict',           # Fail on non-conforming BEAST XML file
	'-beagle_CPU',       # use CPU instance		
	'-beagle_GPU',       # use GPU instance if available
	'-beagle_SSE',       # use SSE extensions if available
	'-beagle_cuda',      # use CUDA parallization if available
	'-beagle_opencl',    # use OpenCL parallization if available
	'-beagle_single',    # use single precision if available
	'-beagle_double',    # use double precision if available	
	'-overwrite',        # Allow overwriting of log files


=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    $self->_set_from_args(
        \@args,
        '-methods' => [ @beast_PARAMS, @beast_SWITCHES ],
        '-create'  => 1,
    );
    my ($out) = $self->SUPER::_rearrange( [qw(OUTFILE_NAME)], @args );
    $self->outfile_name( $out || '' );
    return $self;
}

sub _setparams {
    my ( $self, $file ) = @_;
    my $param_string = '';

	# iterate over parameters and switches
    for my $attr (@beast_PARAMS) {
        my $value = $self->$attr();
        next unless defined $value;
        $param_string .= ' -' . $attr . ' ' . $value;
    }
    for my $attr (@beast_SWITCHES) {
        my $value = $self->$attr();
        next unless $value;
        $param_string .= ' -' . $attr;
    }

    # Set default output file if no explicit output file has been given
    if ( ! $self->outfile_name ) {
        my ( $tfh, $outfile ) = $self->io->tempfile( '-dir' => $self->work_dir );
        close $tfh;
        undef $tfh;
        $self->outfile_name($outfile);
    }
    
    # set input file name
	$param_string .= ' ' . $file;
    
    # hide stderr
    my $null = File::Spec->devnull;
    $param_string .= " > $null 2> $null" if $self->quiet() || $self->verbose < 0;

    return $param_string;
}

=item logfile_name

Getter/setter for the location of a log file that can be loaded in tracer
to assess convergence and ESS

=cut

sub logfile_name {
	my $self = shift;
	$self->{'_logfile_name'} = shift if @_;
	return $self->{'_logfile_name'};
}

=item beastfile_name

Getter/setter for the location of the xml file that
is generated as input for *BEAST

=cut

sub beastfile_name {
	my $self = shift;
	$self->{'_beastfile_name'} = shift if @_;			
	
	# set filename to temporary one if not present
	if ( ! $self->{'_beastfile_name'} ) {
		my ($fh, $fn) = tempfile();
		$self->{'_beastfile_name'} = $fn; 
	}	
	return $self->{'_beastfile_name'};
}

=item log_freq

Getter/setter for the logging frequency.

=cut


sub log_freq {
	my $self = shift;
	$self->{'_log_freq'} = shift if @_;
	return $self->{'_log_freq'};

}

=item chain_length

Returns the number of steps for the Markov chain for running *BEAST.
If chain length is provided it is set to the given value.

=cut

sub chain_length {
	my $self = shift;
	$self->{'_chain_length'} = shift if @_;
	return $self->{'_chain_length'};
}

=item sample_freq

Returns the number determining how often the Markov Chain is sampled.
If the sample frequence is provided it is set to the given value.

=cut

sub sample_freq {
	my $self = shift;
	$self->{'_sample_freq'} = shift if @_;
	return $self->{'_sample_freq'};
}

=item collapse_species

Getter/Setter for flag if taxa of rank subspecies, varietas and formas of one 
species should be collapsed in the inferred species tree

=cut

sub collapse_species {
	my $self = shift;
	$self->{'_collapse_species'} = shift if @_;
	return $self->{'_collapse_species'};
}

=item root_height

Sets the root height to provided value. If no value is given, the root height is
computed from the distance matrix of the provided alignment.
Returns the root height that was set.

=cut 

sub root_height {
	my ( $self, $aln, $value ) = @_;
	if ( defined $value ) {
		$aln->set_namespaces( $beast_pre => $beast_ns );
		$aln->set_meta_object( "${beast_pre}:rootHeight" => $value );
	}
	else {
		$value = $aln->get_meta_object( "${beast_pre}:rootHeight" );
		if ( not defined $value ) {
		
			# XXX !!!
			return $self->root_height( $aln => 1 );
		
		    # write alignment to tempfile to create AlignI
		    my ( $fh, $name ) = tempfile();
		    my @names;
		    $aln->visit(sub{
		    	my $row  = shift;
		    	my $name = $row->get_name;
		    	my $seq  = $row->get_char;
		    	print $fh '>', $name, "\n", $seq, "\n";
		    	push @names, $name;
		    });
			my $alnin = Bio::AlignIO->new( '-format' => 'fasta', '-file' => $name );
			my $alni = $alnin->next_aln;
			
			# compute distance matrix
			my $stats = Bio::Align::DNAStatistics->new();   
			my $jcmatrix = $stats->distance( 
				'-align'  => $alni, 
				'-method' => 'Jukes-Cantor' 
			);
			my $total;
			my $count;
			for my $i ( 0 .. ( $#names - 1 ) ) {
				for my $j ( ( $i + 1 ) .. $#names ) {
					$total += $jcmatrix->get_entry( $names[$i], $names[$j] );
					$count++;
				}
			}

			return $self->root_height( $aln => ( $total / $count ) );
		}
	}
	return $value;
}

=item program_name

Returns the name of the executable.

=cut

sub program_name { $PROGRAM_NAME }

=item program_dir

Not yet implemented.

=cut

sub program_dir { undef }

=item run

Runs *Beast using the data in the provided NeXML file. 
Returns the output file generated from the *Beast run.

=cut

sub run {
	my ($self,$nexml) = @_;
	
	# parse the input file
	$log->info("going to read nexml file $nexml");	
	$self->_alignment( parse(
		'-format' => 'nexml',
		'-file'   => $nexml,
		'-as_project' => 1,
	) );
	
	# generate BEAST xml
	$log->info("going to generate beast xml");	
	my $twig = $self->_make_beast_xml;

	my $filename = $self->beastfile_name();
	$log->info("writing beast input file to $filename");
	open (my $fh, ">", $filename) or die "Error writing BEAST input file: $!"; 
	$twig->print($fh);
	$log->info("beast xml written to $filename");
	close $fh;
	# create the invocation and run the command	
	my $command = $self->executable . $self->_setparams($filename);
	$log->info("going to execute '$command'");	
    my $status  = system $command;
    
    # fetch the output file
    my $outfile = $self->outfile_name();
    if ( !-e $outfile || -z $outfile ) {
        $self->warn("*BEAST call had status of $status: $? [command $command]\n");
        return undef;
    }
	return $outfile;
}

=item version

Returns the version number of *Beast used.

=cut

sub version {
    my ($self) = @_;
    my $exe;
    return undef unless $exe = $self->executable;
    my $string = `$exe -version 2>&1`;
    $string =~ /BEAST (v\d+\.\d+\.\d+)/;
    return version->parse($1) || undef;
}

sub _alignment {
	my ( $self, $thing, $format ) = @_;
	if ( $thing ) {
		if ( -e $thing ) {
			$self->{'_alignment'} = parse(
				'-format'     => $format,
				'-file'       => $thing,
				'-as_project' => 1,
			);
			$self->_validate;
		}
		elsif ( ref $thing ) {
			$self->{'_alignment'} = $thing;
			$self->_validate;
		}
	}
	return $self->{'_alignment'};
}

sub _validate {
	my $self = shift;
	my $project = $self->{'_alignment'};
	my ($taxa) = @{ $project->get_items(_TAXA_) };
	my %taxon = map { $_->get_id => $_ } @{ $taxa->get_entities };
	my $ntax = $taxa->get_ntax;
	
	# add missing rows to sparse matrices
	MATRIX: for my $matrix ( @{ $project->get_items(_MATRIX_) } ) {
		my @rows = @{ $matrix->get_entities };
		#next MATRIX if $ntax == scalar @rows;
		my $nchar = $matrix->get_nchar;
		my %seen;
		$seen{ $_->get_taxon->get_id }++ for @rows;
		my @missing = grep { ! $seen{$_} } keys %taxon;
		my $to = $matrix->get_type_object;
		for my $tid ( @missing ) {
			$matrix->insert( $fac->create_datum( 
				'-type_object' => $to,
				'-name'        => $taxon{$tid}->get_name,
				'-taxon'       => $taxon{$tid},
				'-char'        => ( 'N' x $nchar ),
			) );
		}
	}
}

sub _escape {
	my $string = shift;
	$string =~ s/ /_/g;
	return $string;
}

sub _elt { XML::Twig::Elt->new(@_) }

sub _parse { XML::Twig::Elt->parse(@_) }

sub _make_taxa_xml {
	my ( $self, %args ) = @_;
	my $taxa = _elt( 'taxa', { 'id' => ( $args{'id'} || 'taxa' ) } );
	for my $name ( @{ $args{'taxa'} } ) {
		my $taxon = _elt( 'taxon', { 'id' => _escape($name) . '_taxon' } );
		$taxon->paste($taxa);
	}
	return $taxa;
}

sub _make_sequence_xml {
	my ( $self, $seq ) = @_;
	my $sequence = _elt( 'sequence', $seq->get_char );
	my $name  = $seq->get_name;
	my $taxon = _elt( 'taxon', { 'idref' => _escape($name) . '_taxon' } );
	$taxon->paste($sequence);
	return $sequence;
}

sub _make_alignment_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->get_name;
	my $alignment = _elt( 'alignment', { 'id' => $id, 'dataType' => 'nucleotide' } );
	for my $seq ( @{ $aln->get_entities } ) {
		my $sequence = $self->_make_sequence_xml($seq);
		$sequence->paste( 'last_child' => $alignment );
	}
	return $alignment;
}

sub _make_patterns_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->get_name;
	my $patterns  = _elt( 'patterns',  { 'id'    => "${id}.patterns", 'from' => 1 } );
	my $alignment = _elt( 'alignment', { 'idref' => $id } );
	$alignment->paste($patterns);
	return $patterns;
}

sub _make_constant_size_xml {
	my $self = shift;
	my $constantSize   = _elt( 'constantSize', { 'id' => 'constant', 'units' => 'substitutions' } );
	my $populationSize = _elt( 'populationSize' );
	my $parameter = _elt(
		'parameter' => {
			'id'    => 'constant.popSize',
			'value' => '0.02',
			'lower' => '0.0',
			'upper' => 'Infinity',
		}
	);
	$parameter->paste($populationSize);
	$populationSize->paste($constantSize);
	return $constantSize;
}

sub _make_coalescent_tree_xml {
	my ( $self, $aln, $rootHeight ) = @_;
	my $id = $aln->get_name;
	my $coalescentTree = _elt( 
		'coalescentTree' => { 
			'id'         => "${id}.startingTree",
			'rootHeight' => $rootHeight,
		}
	);
	_elt( 'constantSize', { 'idref' => 'constant' } )->paste($coalescentTree);
	_elt( 'taxa', { 'idref' => 'taxa' } )->paste($coalescentTree);
	return $coalescentTree;
}

sub _make_tree_model_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->get_name;
	my $model = "${id}.treeModel";
	my $treeModel = _elt( 'treeModel', { 'id' => $model } );
	
	# will re-use this
	my ( $param, $height );
	
	# set all internal nodes
	$param = _elt('parameter',{'id'=>"${model}.allInternalNodeHeights"});
	$height = _elt('nodeHeights',{'internalNodes'=>'true','rootNode'=>'true'});
	$param->paste($height);
	$height->paste($treeModel);
	
	# set internal nodes sans root
	$param = _elt('parameter',{'id'=>"${model}.internalNodeHeights"});
	$height = _elt('nodeHeights',{'internalNodes'=>'true'});
	$param->paste($height);
	$height->paste($treeModel);
	
	# set root
	$param = _elt('parameter',{'id'=>"${model}.rootHeight"});
	$height = _elt('rootHeight');
	$param->paste($height);
	$height->paste($treeModel);
	
	# set tree
	_elt('coalescentTree',{'idref'=>"${id}.startingTree"})->paste($treeModel);
	
	return $treeModel;
}

sub _make_strict_clock_branch_rates_xml {
	my ( $self, $aln, %args ) = @_;
	my $id = $aln->get_name;
	my $scbr = _elt( 'strictClockBranchRates', { 'id' => "${id}.branchRates" } );
	my $parm = _elt( 'parameter', { 'id' => "${id}.clock.rate", 'value' => '1.0', %args } );
	my $rate = _elt( 'rate' );
	$parm->paste($rate);
	$rate->paste($scbr);
	return $scbr;
}

sub _make_hky_model_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->get_name;
	my $HKYModel = _elt( 'HKYModel', { 'id' => "${id}.hky" } );
	
	# make frequencies
	my $freq1 = _elt( 'frequencies' );
	$freq1->paste( 'last_child' => $HKYModel );
	my $frequencyModel = _elt( 'frequencyModel', { 'dataType' => 'nucleotide' } );
	$frequencyModel->paste($freq1);
	my $freq2 = _elt( 'frequencies' );
	$freq2->paste( 'last_child' => $frequencyModel );
	my $param1 = _elt( 'parameter', { 
		'id'    => "${id}.frequencies", 
		'value' => '0.25 0.25 0.25 0.25' } );
	$param1->paste( 'last_child' => $freq2 );
	
	# make kappa
	my $kappa = _elt( 'kappa' );
	my $param2 = _elt( 'parameter', {
		'id'    => "${id}.kappa", 
		'value' => '1.0', 
		'lower' => '1.0E-8', 
		'upper' => 'Infinity' } );
	$param2->paste( 'last_child' => $kappa );
	$kappa->paste( 'last_child' => $HKYModel );
	
	# done
	return $HKYModel;
}

sub _make_site_model_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->get_name;
	my $HKYModel = _elt( 'HKYModel', { 'idref' => "${id}.hky" } );
	my $substitutionModel = _elt( 'substitutionModel' );
	my $siteModel = _elt( 'siteModel', { 'id' => "${id}.siteModel" } );
	$HKYModel->paste($substitutionModel);
	$substitutionModel->paste($siteModel);
	return $siteModel;
}

sub _make_tree_likelihood_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->get_name;
	my $tl = _elt( 'treeLikelihood',   { 'id' => "${id}.treeLikelihood", 'useAmbiguities' => 'false' } );
	my $patterns  = _elt( 'patterns',  { 'idref' => "${id}.patterns"  } );
	my $treeModel = _elt( 'treeModel', { 'idref' => "${id}.treeModel" } );
	my $siteModel = _elt( 'siteModel', { 'idref' => "${id}.siteModel" } );
	my $scbr = _elt( 'strictClockBranchRates', { 'idref' => "${id}.branchRates" } );
	$scbr->paste($tl);
	$siteModel->paste($tl);
	$treeModel->paste($tl);
	$patterns->paste($tl);
	return $tl;
}

sub _make_species_xml {
	my ( $self, %args ) = @_; # taxa, treeModel
	my $species = _elt( 'species', { 'id' => 'species' } );
	
	# make geneTrees element
	my $geneTrees = _elt( 'geneTrees', { 'id' => 'geneTrees' } );
	for my $tm ( @{ $args{'treeModel'} } ) {
		my $treeModel = _elt( 'treeModel', { 'idref' => $tm->att('id') } );
		$treeModel->paste($geneTrees);
	}
	$geneTrees->paste($species);
	
	# make species to individual mapping
	my %species;
	for my $taxon ( @{ $args{'taxa'} } ) {
		my $escaped = _escape($taxon);
		
		# collapse species if flag is set. This means we look if there is
		#  a possible second '_' character in the taxon name. If yes,
		#  the taxon name will be assigned to the corresponding species 
		if ( $self->collapse_species ) {
			if ( $escaped =~ /^([^_]+_[^_]+)/ ) {
				my $binomial = $1;
				$species{$binomial} = [] if not $species{$binomial};
				push @{ $species{$binomial} }, $escaped;
			}
		}
		else {
			$species{$escaped} = [$escaped];
		}	
	}

	# make sp elements
	for my $binomial ( keys %species ) {
		my $id = $binomial;
		#$id .= '_species' if scalar @{ $species{$binomial} } == 1;
		my $sp = _elt( 'sp', { 'id' => $id } );
		for my $instance ( @{ $species{$binomial} } ) {
			_elt( 'taxon', { 'idref' => $instance . '_taxon' } )->paste($sp);
		}
		$sp->paste($species);
	}
	return $species;
}

sub _make_species_tree_xml {
	my $speciesTree = _elt( 'speciesTree', { 'id' => 'sptree', 'constantRoot' => 'true' } );
	my $species = _elt( 'species', { 'idref' => 'species' } );
	my $sppSplitPopulations = _elt( 'sppSplitPopulations', { 'value' => '0.02' } );
	my $parameter = _elt( 'parameter', { 'id' => 'speciesTree.splitPopSize' } );
	$parameter->paste($sppSplitPopulations);
	$sppSplitPopulations->paste($speciesTree);
	$species->paste($speciesTree);
	return $speciesTree;
}

sub _make_birth_death_model_xml {
	my $birthDeathModel = _elt( 'birthDeathModel', { 
		'id'    => 'birthDeath',
		'units' => 'substitutions',		
	} );
	my $birthMinusDeathRate = _elt( 'birthMinusDeathRate' );
	$birthMinusDeathRate->paste($birthDeathModel);
	my $parm1 = _elt( 'parameter', {
		'id'    => 'species.birthDeath.meanGrowthRate',
		'value' => '1.0',
		'lower' => '0.0',
		'upper' => 'Infinity',
	} );
	$parm1->paste($birthMinusDeathRate);
	my $relativeDeathRate = _elt( 'relativeDeathRate' );
	$relativeDeathRate->paste( 'last_child' => $birthDeathModel );
	my $parm2 = _elt( 'parameter', {
		'id'    => 'species.birthDeath.relativeDeathRate',
		'value' => '0.5',
		'lower' => '0.0',
		'upper' => '1.0',
	} );
	$parm2->paste($relativeDeathRate);
	return $birthDeathModel;
}

sub _make_speciation_likelihood_xml {
	my $speciationLikelihood = _elt( 'speciationLikelihood', { 
		'id' => 'speciation.likelihood', 
	} );
	my $model = _elt( 'model' );
	$model->paste( $speciationLikelihood );
	my $birthDeathModel = _elt( 'birthDeathModel', {
		'idref' => 'birthDeath',
	} );
	$birthDeathModel->paste( $model );
	my $spec1 = _elt( 'speciesTree' );
	my $spec2 = _elt( 'speciesTree', { 'idref' => 'sptree' } );
	$spec2->paste( $spec1 );
	$spec1->paste( 'last_child' => $speciationLikelihood );
	return $speciationLikelihood;
}

sub _make_tmrca_statistic_xml {
	my ( $self, $species ) = @_;
	my $tmrcaStatistic = _elt( 'tmrcaStatistic', {
		'id'   => 'speciesTree.rootHeight',
		'name' => 'speciesTree.rootHeight',		
	});
	my $mrca = _elt('mrca');
	my $taxa = _elt('taxa');
	$taxa->paste($mrca);
	$mrca->paste($tmrcaStatistic);
	_elt('speciesTree',{'idref'=>'sptree'})->paste($tmrcaStatistic);
	for my $sp ( $species->descendants('sp') ) {
		_elt('sp',{'idref'=>$sp->att('id')})->paste($taxa);
	}
	return $tmrcaStatistic;
}

sub _make_species_coalescent_xml {
	my $speciesCoalescent = _elt( 'speciesCoalescent', { 
		'id' => 'species.coalescent',
	} );
	my $species = _elt( 'species', { 'idref' => 'species' } );
	$species->paste($speciesCoalescent);
	my $speciesTree = _elt( 'speciesTree', { 'idref' => 'sptree' } );
	$speciesTree->paste( 'last_child' => $speciesCoalescent );
	return $speciesCoalescent;
}

sub _make_mixed_distribution_likelihood_xml {
	my ( $self, $species ) = @_;
	my @sp = $species->descendants('sp');
	my $ntax = scalar @sp;
	my @values;
	push @values, 1 for 1 .. $ntax;
	push @values, 0 for 1 .. ( $ntax + ( $ntax - 2 ) );
	
	# root element
	my $mixedDistributionLikelihood = _elt( 'mixedDistributionLikelihood', {
		'id' => 'species.popSize',
	} );
	
	# distribution0 element structure
	my $distribution0 = _elt( 'distribution0' );
	$distribution0->paste( 'last_child' => $mixedDistributionLikelihood );
	my $gammaDistributionModel1 = _elt( 'gammaDistributionModel' );
	$gammaDistributionModel1->paste($distribution0);
	my $shape1 = _elt( 'shape' => 2 );
	$shape1->paste( $gammaDistributionModel1 );
	my $scale1 = _elt( 'scale' );
	$scale1->paste( $gammaDistributionModel1 );
	my $parm1 = _elt( 'parameter', {
		'id'    => 'species.popMean',
		'value' => '0.02',	
	} );
	$parm1->paste($scale1);
	
	# distribution1 element structure
	my $distribution1 = _elt( 'distribution1' );
	$distribution1->paste( 'last_child' => $mixedDistributionLikelihood );
	my $gammaDistributionModel2 = _elt( 'gammaDistributionModel' );
	$gammaDistributionModel2->paste( $distribution1 );
	my $shape2 = _elt( 'shape' => 4 );
	$shape2->paste( $gammaDistributionModel2 );
	my $scale2 = _elt( 'scale' );
	$scale2->paste( 'last_child' => $gammaDistributionModel2 );
	my $parm2 = _elt( 'parameter', {
		'idref' => 'species.popMean',
	} );
	$parm2->paste( $scale2 );
	
	# data element structure
	my $data = _elt( 'data' );
	$data->paste( 'last_child' => $mixedDistributionLikelihood );
	my $parm3 = _elt( 'parameter', {
		'idref' => 'speciesTree.splitPopSize',
	} );
	$parm3->paste( $data );
	
	# indicators element structure
	my $indicators = _elt( 'indicators' );
	$indicators->paste( 'last_child' => $mixedDistributionLikelihood );
	my $parm4 = _elt( 'parameter', { 'value' => "@values" } );
	$parm4->paste( $indicators );
	
	# done
	return $mixedDistributionLikelihood;
}

sub _make_scale_operator_xml {
	my ( $self, %args ) = @_;
	my $factor = $args{'factor'};
	my $weight = $args{'weight'};
	my $param  = $args{'param'};
	my $scaleOperator = _elt( 'scaleOperator', {
		'scaleFactor' => $factor,
		'weight'      => $weight,
	} );
	my $parameter = _elt( 'parameter', { 'idref' => $param } );
	$parameter->paste( $scaleOperator );
	return $scaleOperator;
}

sub _make_delta_exchange_xml {
	my ( $self, %args ) = @_;
	my $delta  = $args{'delta'};
	my $weight = $args{'weight'};
	my $param  = $args{'param'};
	my $deltaExchange = _elt( 'deltaExchange', {
		'delta'  => $delta,
		'weight' => $weight,
	} );
	my $parameter = _elt( 'parameter', { 'idref' => $param } );
	$parameter->paste( $deltaExchange );
	return $deltaExchange;
}

sub _make_treemodel_operators_xml {
	my ( $self, $size, $treeModel ) = @_;
	my @result;
	
	# subtreeSlide
	push @result, _elt('subtreeSlide',{'size'=>$size,'gaussian'=>'true','weight'=>'15'});
	_elt( 'treeModel', { 'idref' => $treeModel } )->paste( $result[-1] );
	
	# narrowExchange
	push @result, _elt( 'narrowExchange', { 'weight' => '15' } );
	_elt( 'treeModel', { 'idref' => $treeModel } )->paste( $result[-1] );

	# wideExchange
	push @result, _elt( 'wideExchange', { 'weight' => '3' } );
	_elt( 'treeModel', { 'idref' => $treeModel } )->paste( $result[-1] );
	
	# wilsonBalding
	push @result, _elt( 'wilsonBalding', { 'weight' => '3' } );
	_elt( 'treeModel', { 'idref' => $treeModel } )->paste( $result[-1] );
	
	# scaleOperator
	push @result, _elt( 'scaleOperator', { 'scaleFactor' => '0.75', 'weight' => '3' } );
	_elt( 'parameter', { 'idref' => "${treeModel}.rootHeight" } )->paste( $result[-1] );

	# uniformOperator
	push @result, _elt( 'uniformOperator', { 'weight' => '30' } );
	_elt('parameter',{'idref'=>"${treeModel}.internalNodeHeights"})->paste($result[-1]);

	return @result;
}

sub _make_updown_operator_xml {
	my ( $self, %args ) = @_;
	my @up   = $args{'up'}   ? @{ $args{'up'} }   : ();
	my @down = $args{'down'} ? @{ $args{'down'} } : ();
	my $factor = $args{'factor'};
	my $weight = $args{'weight'};
	my $speciesTree = $args{'speciesTree'};
	
	# make containing upDownOperator element
	my $upDownOperator = _elt('upDownOperator',{'scaleFactor'=>$factor,'weight'=>$weight});
	
	# create down parameter references
	my $down = _elt('down');
	for my $idref ( @down ) {
		_elt('parameter',{'idref'=>$idref})->paste($down);
	}
	if ( $speciesTree ) {
		_elt('speciesTree',{'idref'=>$speciesTree})->paste($down);
	}
	$down->paste($upDownOperator);
		
	# create up parameter references
	my $up = _elt('up');
	for my $idref ( @up ) {
		_elt('parameter',{'idref'=>$idref})->paste($up);
	}
	$up->paste($upDownOperator);
	return $upDownOperator;
}

sub _make_node_re_height_xml {
	my ( $self, $weight ) = @_;
	my $nodeReHeight = _elt( 'nodeReHeight', { 'weight' => $weight } );
	my $species = _elt( 'species', { 'idref' => 'species' } );
	$species->paste($nodeReHeight);
	my $speciesTree = _elt( 'speciesTree', { 'idref' => 'sptree' } );
	$speciesTree->paste( 'last_child' => $nodeReHeight );
	return $nodeReHeight;
}

sub _make_one_on_x_prior_xml {
	my ( $self, $idref ) = @_;
	my $oneOnXPrior = _elt('oneOnXPrior');
	my $parameter = _elt( 'parameter', { 'idref' => $idref } );
	$parameter->paste($oneOnXPrior);
	return $oneOnXPrior;
}

sub _make_gamma_prior_xml {
	my ( $self, $idref ) = @_;
	my $gammaPrior = _elt('gammaPrior',{
		'shape'  => '0.1',
		'scale'  => '10.0',
		'offset' => '0.0',
	});
	_elt('parameter',{'idref'=>$idref})->paste($gammaPrior);
	return $gammaPrior;
}

sub _make_operators_xml {
	my $self = shift;
	my $operators = _elt( 'operators', { 'id' => 'operators' } );
	my %height = map { $_->get_name => $self->root_height($_) } @{ $self->_alignment->get_matrices };
	my @id = keys %height;

	# scaleOperator
	for my $id ( @id ) {
		$self->_make_scale_operator_xml(
			'factor' => '0.75',
			'weight' => '1',
			'param'  => "${id}.kappa"
		)->paste( 'last_child' => $operators );
		$self->_make_delta_exchange_xml(
			'delta'  => '0.01',
			'weight' => 1,
			'param'  => "${id}.frequencies",
		)->paste( 'last_child' => $operators );
	}
		
	# scaleOperator
	for my $id ( @id ) {
		$self->_make_scale_operator_xml(
			'factor' => '0.75',
			'weight' => '3',
			'param'  => "${id}.clock.rate",
		)->paste( 'last_child' => $operators );
	}	
	
	# upDownOperator
	my @up = map { "${_}.clock.rate" } @id;
	push @up, 'species.birthDeath.meanGrowthRate';
	my @down = qw(species.popMean speciesTree.splitPopSize);
	push @down, "${_}.treeModel.allInternalNodeHeights" for @id;
	$self->_make_updown_operator_xml(
		'up'          => [ reverse @up   ],
		'down'        => [ reverse @down ],
		'factor'      => '0.75',
		'weight'      => '30',
		'speciesTree' => 'sptree',
	)->paste( 'last_child' => $operators );
	
	# tree operators
	for my $id ( @id ) {
		my $tm = "${id}.treeModel";
		for my $elt ( $self->_make_treemodel_operators_xml( $height{$id} / 10 , $tm ) ) {
			$elt->paste( 'last_child' => $operators );
		}
	}
	
	# upDownOperator
	my $i = 0;
	for my $id ( @id ) {
		my @up;
		@up = ( "${id}.clock.rate" ) if $i;
		$self->_make_updown_operator_xml(
			'factor' => '0.75',
			'weight' => '3',
			'down'   => [ "${id}.treeModel.allInternalNodeHeights" ],
			'up'     => \@up,
		)->paste( 'last_child' => $operators );
		$i++;
	}
		
	# scaleOperator species.popMean
	$self->_make_scale_operator_xml(
		'factor' => '0.9',
		'weight' => '5',
		'param'  => 'species.popMean',
	)->paste( 'last_child' => $operators );		
				
	# scaleOperator species.birthDeath.meanGrowthRate
	$self->_make_scale_operator_xml(
		'factor' => '0.75',
		'weight' => '3',
		'param'  => 'species.birthDeath.meanGrowthRate',
	)->paste( 'last_child' => $operators );
	
	# scaleOperator species.birthDeath.relativeDeathRate
	$self->_make_scale_operator_xml(
		'factor' => '0.75',
		'weight' => '3',
		'param'  => 'species.birthDeath.relativeDeathRate',
	)->paste( 'last_child' => $operators );
	
	# scaleOperator speciesTree.splitPopSize
	$self->_make_scale_operator_xml(
		'factor' => '0.5',
		'weight' => '94',
		'param'  => 'speciesTree.splitPopSize',
	)->paste( 'last_child' => $operators );
	
	# nodeReHeight
	$self->_make_node_re_height_xml('94')->paste( 'last_child' => $operators );
		
	return $operators;
}

sub _make_prior_xml {
	my $self = shift;
	my $prior = _elt( 'prior', { 'id' => 'prior' } );
	
	# speciesCoalescent "species.coalescent"
	_elt('speciesCoalescent',{'idref'=>'species.coalescent'})->paste('last_child'=>$prior);
	
	# mixedDistributionLikelihood "species.popSize"
	_elt('mixedDistributionLikelihood',{'idref'=>'species.popSize'})->paste('last_child'=>$prior);
	
	# speciationLikelihood "speciation.likelihood"
	_elt('speciationLikelihood',{'idref'=>'speciation.likelihood'})->paste('last_child'=>$prior);
	
	# oneOnXPrior for kappa
	my @id = map { $_->get_name } @{ $self->_alignment->get_matrices };
	$self->_make_one_on_x_prior_xml("${_}.kappa")->paste('last_child'=>$prior) for @id;
	
	# gammaPrior for clock.rate
	$self->_make_gamma_prior_xml("${_}.clock.rate")->paste('last_child'=>$prior) for @id;
	
	# oneOnXPrior for species.popMean
	$self->_make_one_on_x_prior_xml('species.popMean')->paste('last_child'=>$prior);
	
	# oneOnXPrior for species.birthDeath.meanGrowthRate
	$self->_make_one_on_x_prior_xml('species.birthDeath.meanGrowthRate')->paste('last_child'=>$prior);
	
	return $prior;

}

sub _make_likelihood_xml {
	my $self = shift;
	my $likelihood = _elt( 'likelihood', { 'id' => 'likelihood' } );
	my @id = map { $_->get_name } @{ $self->_alignment->get_matrices };
	my $tl = 'treeLikelihood';
	for my $id ( @id ) {
		_elt( $tl, { 'idref' => "${id}.${tl}" } )->paste( 'last_child' => $likelihood );
	}
	return $likelihood;
}

sub _make_posterior_xml {
	my $self = shift;
	my $posterior = _elt( 'posterior', { 'id' => 'posterior' } );
	$self->_make_prior_xml->paste(      'last_child' => $posterior );
	$self->_make_likelihood_xml->paste( 'last_child' => $posterior );
	return $posterior;
}

sub _make_species_logtree_xml {
	my $self = shift;
	my $logTree = _elt( 'logTree', {
		'id'                   => 'species.treeFileLog',
		'logEvery'             => $self->sample_freq,
		'nexusFormat'          => 'true',
		'fileName'             => $self->outfile_name,
		'sortTranslationTable' => 'true',
	} );
	_elt( 'speciesTree', { 'idref' => 'sptree' }  )->paste( 'last_child' => $logTree );
	_elt( 'posterior', { 'idref' => 'posterior' } )->paste( 'last_child' => $logTree );
	return $logTree;
}

sub _make_logfile_xml {
	my $self = shift;
	my $file = $self->logfile_name;
	my $freq = $self->log_freq || 1000;
	
	# create enclosing element
	my $log = _elt( 'log', {
		'logEvery' => $freq,
		'fileName' => $file,
	} );
	
	# create static columns
	my %cols = (
		'posterior'            => 'posterior',
		'prior'                => 'prior',
		'likelihood'           => 'likelihood',
		'speciesCoalescent'    => 'species.coalescent',
		'speciationLikelihood' => 'speciation.likelihood',
		'tmrcaStatistic'       => 'speciesTree.rootHeight',
	);
	for my $col ( keys %cols ) {
		_elt( $col, { 'idref' => $cols{$col} } )->paste( 'last_child' => $log );
	}
	
	# create static parameters
	my @params = qw(
		species.popMean 
		speciesTree.splitPopSize 
		species.birthDeath.meanGrowthRate 
		species.birthDeath.relativeDeathRate
	);
	for my $param ( @params ) {
		_elt( 'parameter', { 'idref' => $param } )->paste( 'last_child' => $log );
	}
	
	# create dynamic parameters
	my @dynparam = qw(
		treeModel.rootHeight
		kappa
		frequencies
		clock.rate
	);
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		my $id = $aln->get_name;
		for my $param ( @dynparam ) {
			my $idref = "${id}.${param}";
			_elt( 'parameter', { 'idref' => $idref } )->paste( 'last_child' => $log );
		}
	}
	
	# create tree likelihoods
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		my $id = $aln->get_name;
		my $tl = 'treeLikelihood';
		my $idref = "${id}.${tl}";
		_elt( $tl, { 'idref' => $idref } )->paste( 'last_child' => $log );
	}	
	
	# done
	return $log;
}

sub _make_mcmc_xml {
	my $self = shift;
	my $mcmc = _elt( 'mcmc', {
		'id'           => 'mcmc', 
		'chainLength'  => $self->chain_length, 
		'autoOptimize' => 'true',
	} );
	$self->_make_posterior_xml->paste( 'last_child' => $mcmc );
	_elt( 'operators', { 'idref' => 'operators' } )->paste( 'last_child' => $mcmc );
	if ( $self->logfile_name ) {
		$self->_make_logfile_xml->paste( 'last_child' => $mcmc );
	}	
	$self->_make_species_logtree_xml->paste( 'last_child' => $mcmc );
	return $mcmc;
}

sub _make_beast_xml {
	my $self = shift;
	my $twig = XML::Twig->new( 'pretty_print' => 'indented' );
	my $root = _elt('beast');
	$twig->set_root($root);
	my %paste = ( 'last_child' => $root );
	
	# gather binomial names across all alignments, make taxa element
	my %labels;
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		for my $seq ( @{ $aln->get_entities } ) {
			my $name = $seq->get_name;
			$labels{$name}++;
		}
	}
	my @labels = keys %labels;	
	$self->_make_taxa_xml('taxa'=>\@labels)->paste(%paste);
	
	# make alignment elements
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		$self->_make_alignment_xml($aln)->paste(%paste);
	}
	
	# make pattern elements
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		$self->_make_patterns_xml($aln)->paste(%paste);
	}
	
	# make constantSize element
	$self->_make_constant_size_xml->paste(%paste);
	
	# make coalescentTree elements
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		my $height = $self->root_height($aln);
		$self->_make_coalescent_tree_xml( $aln => $height )->paste(%paste);
	}
	
	# make treeModel elements
	my @treeModel;
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		push @treeModel, $self->_make_tree_model_xml($aln);
		$treeModel[-1]->paste(%paste);
	}
	
	# make strictClockBranchRates elements
	my $i = 0;
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		my %args;
		%args = ( 'lower' => '0.0', 'upper' => 'Infinity' ) if $i;
		$self->_make_strict_clock_branch_rates_xml($aln, %args)->paste(%paste);
		$i++;
	}
	
	# make HKYModel elements
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		$self->_make_hky_model_xml($aln)->paste(%paste);
	}
	
	# make siteModel elements
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		$self->_make_site_model_xml($aln)->paste(%paste);
	}
	
	# make treeLikelihood elements
	for my $aln ( @{ $self->_alignment->get_matrices } ) {
		$self->_make_tree_likelihood_xml($aln)->paste(%paste);
	}
	
	# make species element
	my $species = $self->_make_species_xml('taxa'=>\@labels,'treeModel'=>\@treeModel);
	$species->paste(%paste);
	
	# make speciesTree element
	$self->_make_species_tree_xml->paste(%paste);
	
	# make birthDeathModel element
	$self->_make_birth_death_model_xml->paste(%paste);
	
	# make speciationLikelihood element
	$self->_make_speciation_likelihood_xml->paste(%paste);
	
	# make tmrcaStatistic element
	$self->_make_tmrca_statistic_xml($species)->paste(%paste);
	
	# make speciesCoalescent element
	$self->_make_species_coalescent_xml->paste(%paste);
	
	# make mixedDistributionLikelihood element
	$self->_make_mixed_distribution_likelihood_xml($species)->paste(%paste);
	
	# make operators element
	$self->_make_operators_xml->paste(%paste);
	
	# make mcmc element
	$self->_make_mcmc_xml->paste(%paste);
	
	return $twig;
}

=back

=cut

1;
