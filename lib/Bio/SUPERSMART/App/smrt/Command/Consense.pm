package Bio::SUPERSMART::App::smrt::Command::Consense;
use strict;
use warnings;
use File::Temp 'tempfile';
use Bio::Phylo::Factory;
use Bio::Phylo::PhyLoTA::Service::TreeService;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

my $ts  = Bio::Phylo::PhyLoTA::Service::TreeService->new;
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
	my $intree_default = "chronogram.dnd";
	return (
		["infile|i=s", "newick input tree(s) file", { arg => "file", default => $intree_default}],
		["outfile|o=s", "name of the output file, defaults to '$outfile_default'", {default => $outfile_default, arg => "file"}],
		["format|f=s", "format of consensus tree file, (nexus, newick) defaults to 'nexus'", {default => 'nexus'}],
        ["burnin|b=f", "fraction of burnin to omit", {default => $ts->config->BURNIN}],
        ["heights|h=s", "how to summarize heights (keep, median, mean, ca)", {default => $ts->config->NODE_HEIGHTS}],
        ["limit|l=f", "the minimum support for a node to be annotated",{default => 0.0}],
        ["prob|p","write node support as probabilities (otherwise fractions)",{}],
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
	$self->logger->info("preparing input file");
    my ( $fh, $filename ) = tempfile();
    close $fh;
    my $ntrees = $ts->newick2nexus( $opt->infile => $filename );

    # run treeannotator, remove temp file
    $self->logger->info("computing consensus");
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
        if ( my $posterior = $node->get_meta_object('fig:posterior') ) {

            # apply optionally converted node support as node label
			if ( $opt->prob ) {
				$node->set_name( $posterior );
			}
			else {
				my $count = int( $posterior * $ntrees + 0.5 );
				$node->set_name( "$count/$ntrees");
			}
		}
    });

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
