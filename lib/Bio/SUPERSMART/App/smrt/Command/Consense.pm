package Bio::SUPERSMART::App::smrt::Command::Consense;
use strict;
use warnings;
use File::Temp 'tempfile';
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'unparse';
use Bio::SUPERSMART::Service::TreeService;
use Bio::SUPERSMART::Service::CalibrationService;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

my $ts  = Bio::SUPERSMART::Service::TreeService->new;
my $fac = Bio::Phylo::Factory->new;

# ABSTRACT: constructs a consensus tree

=head1 NAME

Consense.pm - constructs a consensus tree

=head1 SYNOPSYS

smrt consense [-v] [-w <dir>] -i <file> [-o <file>] [-f <format>] [-b <burnin>] 
[-h <heights>] [-l <limit>] [-p]

=head1 DESCRIPTION

Given an input set of trees in newick format, creates an annotated consensus tree.

=cut


sub options {
	my ($self, $opt, $args) = @_;
	my $outfile_default = "consensus.nex";
	my $intree_default  = "chronogram.dnd";
	my $fossils_default = "fossils.tsv";
	my $heights_default = $ts->config->NODE_HEIGHTS; 
	return (
		[
		     "infile|i=s", 
    		 "newick input tree(s) file", 
		     { arg => "file", default => $intree_default, galaxy_in => 1, galaxy_type => 'data'}
		],
		[
		     "outfile|o=s", 
		     "name of the output file, defaults to '$outfile_default'", 
		     { default => $outfile_default, arg => "file", galaxy_out => 1, galaxy_type => 'data', galaxy_label => 'consensus'}
		],
		[
		     "burnin|b=s", "fraction of burnin to omit, defaults to " . $ts->config->BURNIN . "; set to 0.0 for no burnin", 
    		 { default => $ts->config->BURNIN, arg => "fraction",  galaxy_in => 1, galaxy_type => 'text'}
		],
		[
		     "heights|e=s", 
		     "how to summarize heights (keep, median, mean, ca), defaults to $heights_default", 
		     { default => $heights_default, arg => "keep|median|mean|ca", 
			   galaxy_in => 1, galaxy_type => 'select', galaxy_value => $heights_default, galaxy_options => ['keep', 'median', 'mean', 'ca']}
		],
		[
		     "limit|l=f", 
		     "the minimum support for a node to be annotated",
		     { default => 0.0, arg => "value", galaxy_in => 1, galaxy_type => 'text', galaxy_value => '0.0'}
		],
		[
		     "fossils|s=s",
		     "fossil table (if re-applying calibration points)", 
		     {default => $fossils_default, arg => "file", galaxy_in => 1, galaxy_type => 'data' }
		],
		[
		     "format|f=s", 
		     "format of consensus tree file, (nexus, newick) defaults to 'nexus'", 
		     { default => 'nexus', arg => "format" },
		],
		[    
			 "prob|p",
			 "write node support as probabilities (otherwise fractions)",
			 { galaxy_in => 1, galaxy_type => 'boolean'}
		],
	);
}

sub validate {
	my ($self, $opt, $args) = @_;

	my $file = $opt->infile;
	$self->usage_error("no infile argument given")     if not $file;
	$self->usage_error("file $file does not exist") unless -e $file;
	$self->usage_error("file $file is empty")       unless -s $file;

	$self->usage_error("only newick or nexus for now")     if $opt->format !~ /^(?:newick|nexus)$/i;
	$self->usage_error("burnin should be between 0 and 1") if $opt->burnin < 0 or $opt->burnin > 1;
	$self->usage_error("limit should be between 0 and 1")  if $opt->limit  < 0 or $opt->limit  > 1;
	$self->usage_error("heights should be one of 'keep', 'median', 'mean' or 'ca'") if $opt->heights !~ /(?:keep|median|mean|ca)/;
}

sub run {
	my ($self, $opt, $args) = @_;

	# write simple nexus file (no translation table)
	$self->logger->info("Preparing input file");
	my ( $fh, $filename ) = tempfile();
	close $fh;
	my $ntrees = $ts->newick2nexus( $opt->infile => $filename );

	# run treeannotator, remove temp file
	$self->logger->info("Computing consensus with heights '".$opt->heights."'");
	my $consensus = $ts->consense_trees(
		'-infile'  => $filename,
		'-burnin'  => $opt->burnin,
		'-heights' => $opt->heights,
		'-limit'   => $opt->limit,
	);
	unlink $filename;

	# update node labels to distinguish bootstraps from posteriors
	$self->logger->info("applying node labels");
	$consensus->visit(sub{
        	my $node = shift;
		if ( $node->is_internal ) {
			if ( my $posterior = $node->get_meta_object('fig:posterior') ) {
	
				# apply optionally converted node support as node label
				if ( not $opt->prob ) {
					my $count = int( $posterior * $ntrees + 0.5 );
					$node->set_meta_object( 'fig:bootstrap' => $count );
					my @meta = @{ $node->get_meta('fig:posterior') };
					$node->remove_meta( $_ ) for @meta;
				}
			}
			else {
				$self->logger->debug("no posterior probability on node");
				$self->logger->debug($node->to_js);
			}
		}
	});

	# apply calibration points
	if ( $opt->fossils and -e $opt->fossils ) {
		$ts->remap_to_ti($consensus);
		my $cs = Bio::SUPERSMART::Service::CalibrationService->new;
		my @fossils = $cs->read_fossil_table( $opt->fossils );
		my @identified = map { $cs->find_calibration_point($_) } @fossils;
		my $ct = $cs->create_calibration_table( $consensus, @identified );
		$ts->remap_to_name($consensus);
	}

	# generate output
	my %args = ( '-nodelabels' => 1 );
	my $string;
	if ( $opt->format =~ /nexus/i ) {

		# create and populate a mesquite-like nexus project
		my $project = $fac->create_project;
		my $forest  = $fac->create_forest;
		$forest->insert($consensus);
		my $taxa = $forest->make_taxa;
		$taxa->set_forest($forest);
		$project->insert($taxa);
		$project->insert($forest);
		$string = unparse(
        		'-format' => 'figtree',
        		'-phylo'  => $project,
		);
	}
	else {
		$string = $consensus->to_newick(%args);
	}
    
	# write to file
	open my $out, '>', $self->outfile or die $!;
	print $out $string;
	$self->logger->info("DONE, results written to ".$self->outfile);
}

1;
