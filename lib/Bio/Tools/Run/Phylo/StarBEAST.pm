package Bio::Tools::Run::Phylo::StarBEAST;
use strict;
use XML::Twig;
use base qw(Bio::Root::Root Bio::Tools::Run::PhyloBase);

our $PROGRAM_NAME = 'beast';
our @beast_PARAMS = qw(mc3_chains mc3_delta mc3_temperatures mc3_swap seed threshold);
our @beast_SWITCHES = qw(verbose warnings strict);

sub alignments {
	my $self = shift;
	if ( @_ ) {
		$self->{alignments} = \@_;
	}
	ref $self->{alignments} ? @{ $self->{alignments} } : ();
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
		my $taxon = _elt( 'taxon', { 'id' => _escape($name) } );
		$taxon->paste($taxa);
	}
	return $taxa;
}

sub _make_sequence_xml {
	my ( $self, $seq ) = @_;
	my $sequence = _elt( 'sequence', $seq->seq );
	my $name  = $seq->species->binomial('FULL');
	my $taxon = _elt( 'taxon', { 'idref' => _escape($name) } );
	$taxon->paste($sequence);
	return $sequence;
}

sub _make_alignment_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->id;
	my $alignment = _elt( 'alignment', { 'id' => $id, 'dataType' => 'nucleotide' } );
	for my $seq ( $aln->each_seq ) {
		my $sequence = $self->_make_sequence_xml($seq);
		$sequence->paste( 'last_child' => $alignment );
	}
	return $alignment;
}

sub _make_patterns_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->id;
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
	my $id = $aln->id;
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
	my $id = $aln->id;
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

sub _make_string_clock_branch_rates_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->id;
	my $scbr = _elt( 'strictClockBranchRates', { 'id' => "${id}.branchRates" } );
	my $parm = _elt( 'parameter', { 'id' => "${id}.clock.rate", 'value' => '1.0' } );
	my $rate = _elt( 'rate' );
	$parm->paste($rate);
	$rate->paste($scbr);
	return $scbr;
}

sub _make_hky_model_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->id;
_parse(<<"HERE"
<HKYModel id="${id}.hky">
	<frequencies>
		<frequencyModel dataType="nucleotide">
			<frequencies>
				<parameter id="${id}.frequencies" value="0.25 0.25 0.25 0.25"/>
			</frequencies>
		</frequencyModel>
	</frequencies>
	<kappa>
		<parameter id="${id}.kappa" value="1.0" lower="1.0E-8" upper="Infinity"/>
	</kappa>
</HKYModel>
HERE
)
}

sub _make_site_model_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->id;
	my $HKYModel = _elt( 'HKYModel', { 'idref' => "${id}.hky" } );
	my $substitutionModel = _elt( 'substitutionModel' );
	my $siteModel = elt( 'siteModel', { 'id' => "${id}.siteModel" } );
	$HKYModel->paste($substitutionModel);
	$substitutionModel->paste($siteModel);
	return $siteModel;
}

sub _make_tree_likelihood_xml {
	my ( $self, $aln ) = @_;
	my $id = $aln->id;
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
		if ( $escaped =~ /^([^_]+_[^_]+)/ ) {
			my $binomial = $1;
			$species{$binomial} = [] if not $species{$binomial};
			push @{ $species{$binomial} }, $escaped;
		}
	}

	# make sp elements
	for my $binomial ( keys %species ) {
		my $id = $binomial;
		$id .= '_species' if scalar @{ $species{$binomial} } == 1;
		my $sp = _elt( 'sp', { 'id' => $id } );
		for my $instance ( @{ $species{$binomial} } ) {
			_elt( 'taxon', { 'idref' => $instance } )->paste($sp);
		}
		$sp->paste($species);
	}
	return $species;
}

sub _make_species_tree_xml {
_parse(<<"HERE"
	<speciesTree id="sptree" constantRoot="true">
		<species idref="species"/>
		<sppSplitPopulations value="0.02">
			<parameter id="speciesTree.splitPopSize"/>
		</sppSplitPopulations>
	</speciesTree>
HERE
)
}

sub _make_birth_death_model_xml {
_parse(<<'HERE'
<birthDeathModel id="birthDeath" units="substitutions">
	<birthMinusDeathRate>
		<parameter 
			id="species.birthDeath.meanGrowthRate" 
			value="1.0" 
			lower="0.0" 
			upper="Infinity"/>
	</birthMinusDeathRate>
	<relativeDeathRate>
		<parameter 
			id="species.birthDeath.relativeDeathRate" 
			value="0.5" 
			lower="0.0" 
			upper="1.0"/>
	</relativeDeathRate>
</birthDeathModel>
HERE
)
}

sub _make_speciation_likelihood_xml {
_parse(<<'HERE'
<speciationLikelihood id="speciation.likelihood">
	<model>
		<birthDeathModel idref="birthDeath"/>
	</model>
	<speciesTree>
		<speciesTree idref="sptree"/>
	</speciesTree>
</speciationLikelihood>
HERE
)
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
_parse(<<'HERE'
<speciesCoalescent id="species.coalescent">
	<species idref="species"/>
	<speciesTree idref="sptree"/>
</speciesCoalescent>
HERE
)
}

sub _make_mixed_distribution_likelihood_xml {
	my ( $self, $species ) = @_;
	my @sp = $species->descendants('sp');
	my $ntax = scalar @sp;
	my @values;
	push @values, 1 for 1 .. $ntax;
	push @values, 0 for 1 .. ( $ntax + ( $ntax - 2 ) );
_parse(<<"HERE"
	<mixedDistributionLikelihood id="species.popSize">
		<distribution0>
			<gammaDistributionModel>
				<shape>
					2
				</shape>
				<scale>
					<parameter id="species.popMean" value="0.02"/>
				</scale>
			</gammaDistributionModel>
		</distribution0>
		<distribution1>
			<gammaDistributionModel>
				<shape>
					4
				</shape>
				<scale>
					<parameter idref="species.popMean"/>
				</scale>
			</gammaDistributionModel>
		</distribution1>
		<data>
			<parameter idref="speciesTree.splitPopSize"/>
		</data>
		<indicators>
			<parameter value="@values"/>
		</indicators>
	</mixedDistributionLikelihood>
HERE
);
}

sub _make_scale_operator_xml {
	my ( $self, %args ) = @_;
	my $factor = $args{'factor'};
	my $weight = $args{'weight'};
	my $param  = $args{'param'};
	_parse(qq(
		<scaleOperator scaleFactor="$factor" weight="$weight">
			<parameter idref="$param"/>
		</scaleOperator>
	));
}

sub _make_delta_exchange_xml {
	my ( $self, %args ) = @_;
	my $delta  = $args{'delta'};
	my $weight = $args{'weight'};
	my $param  = $args{'param'};
	return _parse(qq(
		<deltaExchange delta="$delta" weight="$weight">
			<parameter idref="$param"/>
		</deltaExchange>
	));
}

sub _make_treemodel_operators_xml {
	my ( $self, $size, $treeModel ) = @_;
	return 
		_parse(qq(
		<subtreeSlide size="$size" gaussian="true" weight="15">
			<treeModel idref="$treeModel"/>
		</subtreeSlide>)),
		_parse(qq(
		<narrowExchange weight="15">
			<treeModel idref="$treeModel"/>
		</narrowExchange>)),
		_parse(qq(
		<wideExchange weight="3">
			<treeModel idref="$treeModel"/>
		</wideExchange>)),
		_parse(qq(
		<wilsonBalding weight="3">
			<treeModel idref="$treeModel"/>
		</wilsonBalding>)),
		_parse(qq(
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="${treeModel}.rootHeight"/>
		</scaleOperator>)),
		_parse(qq(		
		<uniformOperator weight="30">
			<parameter idref="${treeModel}.internalNodeHeights"/>
		</uniformOperator>));
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
	_parse(qq(
		<nodeReHeight weight="$weight">
			<species idref="species"/>
			<speciesTree idref="sptree"/>
		</nodeReHeight>));
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

sub _make_beast_xml {
	my $self = shift;
	my $twig = XML::Twig->new( 'pretty_print' => 'indented' );
	my $root = _elt('beast');
	$twig->set_root($root);
	
	# gather binomial names across all alignments, make taxa element
	my @labels;
	for my $aln ( $self->alignments ) {
		for my $seq ( $aln->each_seq ) {
			push @labels, $seq->species->binomial('FULL');
		}
	}
	$self->_make_taxa_xml('taxa'=>\@labels)->paste($root);
	
	# make alignment elements
	my @alignments;
	for my $aln ( $self->alignments ) {
		push @alignments, $self->_make_alignment_xml($aln);
	}
	
	# make pattern elements
	my @patterns;
	for my $aln ( $self->alignments ) {
		push @patterns, $self->_make_patterns_xml($aln);
	}
	
	# make constantSize element
	my $constantSize = $self->_make_constant_size_xml;
	
	# make coalescentTree elements
	my @coalescentTree;
	for my $aln ( $self->alignments ) {
		push @coalescentTree, $self->_make_coalescent_tree_xml($aln);
	}
	
	# make treeModel elements
	my @treeModel;
	for my $aln ( $self->alignments ) {
		push @treeModel, $self->_make_tree_model_xml($aln);
	}
	
	# make strictClockBranchRates elements
	my @strictClockBranchRates;
	for my $aln ( $self->alignments ) {
		push @strictClockBranchRates, $self->_make_string_clock_branch_rates_xml($aln);
	}
	
	# make HKYModel elements
	my @HKYModel;
	for my $aln ( $self->alignments ) {
		push @HKYModel, $self->_make_hky_model_xml($aln);
	}
	
	# make siteModel elements
	my @siteModel;
	for my $aln ( $self->alignments ) {
		push @siteModel, $self->_make_site_model_xml($aln);
	}
	
	# make treeLikelihood elements
	my @treeLikelihood;
	for my $aln ( $self->alignments ) {
		push @treeLikelihood, $self->_make_tree_likelihood_xml($aln);
	}
	
	# make species element
	my $species = $self->_make_species_xml('taxa'=>\@labels,'treeModel'=>\@treeModel);
	
	# make speciesTree element
	my $speciesTree = $self->_make_species_tree_xml;
	
	# make birthDeathModel element
	my $birthDeathModel = $self->_make_birth_death_model_xml;
	
	# make speciationLikelihood element
	my $speciationLikelihood = $self->_make_speciation_likelihood_xml;
	
	# make tmrcaStatistic element
	my $tmrcaStatistic = $self->_make_tmrca_statistic_xml;
	
	# make speciesCoalescent element
	my $speciesCoalescent = $self->_make_species_coalescent_xml;
	
	# make mixedDistributionLikelihood element
	my $mixedDistributionLikelihood = $self->_make_mixed_distribution_likelihood_xml;
}

1;