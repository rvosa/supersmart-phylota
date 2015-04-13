package Bio::SUPERSMART::App::smrt::Command::Pipeline;

use strict;
use warnings;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: Run the entire SUPERSMART pipeline

=head1 NAME

Pipeline.pm - This command will run the whole SUPERSMART pipeline on given input data

=head1 SYNOPSYS

smrt  taxize [-v ] -t <taxa> [-w <dir>] [-f <file>]

=head1 DESCRIPTION

Given one or more higher root taxa, this command will run all steps of the SUPERSMART pipeline in sequence.
Multiple taxa can be entered by separating them with commata.
Optionally, a fossil table can be given as argument and the tree will be calibrated. 
If no working directory is specified via the command-line, this subcommand will create a working directory
with directory name <DATE>-<root_taxon>, where root_taxon is the first one if multiple are given as input. 

Given the complexity of the individual pipeline steps,one should be very careful to run the 
whole pipeline as a 'black box'. However, all intermediate files are stored in the working directory and a logfile
for each individual subcommand is created, so the user can check all intermediate results. 

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	return(
		[
			"root_taxa|t=s",
			"one or multiple taxon names (seperated by commata) to be expanded to rank 'species'. Taxon names containing spaces must be enclosed in quotes",
			{ mandatory=>1 }
		],

		[
			"outgroup|o=s",
			"one or multiple taxon names (seperated by commata) to be included of outgroup. Taxon names containing spaces must be enclosed in quotes",
		],

		[
			"fossiltable|f=s", 
			"tsv (tab-separated value) file containing fossil table with at least 5 columns (id, name, crown/stem, taxon, age)", 
			{ arg => "file"}
		]			
	)
}

sub validate {
	my ($self, $opt, $args) = @_;		
	if ( ! $opt->root_taxa ){
		$self->usage_error("need root_taxa argument");
	}
}


sub run{
	my ($self, $opt, $args) = @_;		
	my $log = $self->logger;

	my $root_taxa = $opt->root_taxa;
	my $outgroup = $opt->outgroup;
	$root_taxa .= "," . $outgroup if $outgroup;
	my $workdir = $opt->workdir;	
	
	# if no workdir is given, create directory with name being a combination of date and first root taxon	
	if ( ! $workdir  ) {		
		my ($taxon) = split(',', $root_taxa);				
		require POSIX; 
		require Cwd;
		my $date = POSIX::strftime "%m-%d-%Y", localtime;		
		$workdir = Cwd::getcwd . '/' . $date . '-' . $taxon;
		$self->workdir($workdir);
		$log->info("No working directory given in option. Creating working directory $workdir");
		system("mkdir $workdir") and $log->warn($?);
	}
	
	# define file names for intermediate files
	my $taxafile = $self->absolute_path("taxa.tsv");
	my $classtreefile = $self->absolute_path("classtree.dnd");
	my $alnfile = $self->absolute_path("aligned.txt");
	my $orthfile = $self->absolute_path("clusters.txt");
	my $supermatrixfile = $self->absolute_path("supermatrix.phy");
	my $backbonefile = $self->absolute_path("backbone.dnd");
	my $rerootedfile = $self->absolute_path("backbone-rerooted.dnd");
	my $chronogramfile = $self->absolute_path("chronogram.dnd");
	my $finaltreefile = $self->absolute_path("final.dnd");
	
	# other parameters
	my $beast_gens = 100000;
	
	# Step 1: Taxize 
	$log->info("Running pipeline step #1: smrt taxize");
	my $success = $self->_run_subcommand( ("taxize",  "-w", $workdir, "-r", $root_taxa, "-o", $taxafile, "-e", "Species", "-w", $workdir, "-l", "taxize.log") );	
	die ("Pipeline step 'taxize' failed. Try to run 'smrt taxize' manually with '-v' option to find out what went wrong.") unless ( $success and  -e $taxafile and -s $taxafile);
	$log->info("Step #1 smrt taxize succeed");
	
	# Step 2: Classify
	$log->info("Running pipeline step #2: smrt classify");
	$success = $self->_run_subcommand( ("classify",  "-w", $workdir, "-i", $taxafile, "-o", $classtreefile, "-l", "classify.log") );	
	die ("Pipeline step 'classify' failed. Try to run 'smrt classify' manually with '-v' option to find out what went wrong.") unless ( $success and  -e $classtreefile and -s $classtreefile);
	$log->info("Step #2 smrt classify succeed");
	
	# Step 3: Align
	$log->info("Running pipeline step #3: smrt align");
	$success = $self->_run_subcommand( ("align",  "-w", $workdir, "-i", $taxafile, "-o", $alnfile, "-l", "align.log") );	
	die ("Pipeline step 'align' failed. Try to run 'smrt align' manually with '-v' option to find out what went wrong.") unless ( $success and  -e $alnfile and -s $alnfile);
	$log->info("Step #3 smrt align succeed");
	
	# Step 4: Orthologize
	$log->info("Running pipeline step #4: smrt orthologize");
	$success = $self->_run_subcommand( ("orthologize",  "-w", $workdir, "-i", $alnfile, "-o", $orthfile, "-l", "orthologize.log") );	
	die ("Pipeline step 'orthologize' failed. Try to run 'smrt orthologize' manually with '-v' option to find out what went wrong.") unless ( $success and  -e $orthfile and -s $orthfile);
	$log->info("Step #4 smrt orthologize succeed");
	
	# Step 5: BBmerge
	$log->info("Running pipeline step #5: smrt bbmerge");
	$success = $self->_run_subcommand( ("bbmerge",  "-w", $workdir, "-a", $orthfile, "-t", $taxafile, "-o", $supermatrixfile, "-l", "bbmerge.log") );	
	die ("Pipeline step 'bbmerge' failed. Try to run 'smrt bbmerge' manually with '-v' option to find out what went wrong.") unless ( $success and  -e $supermatrixfile and -s $supermatrixfile);
	$log->info("Step #5 smrt bbmerge succeed");
		
	# Step 6: BBinfer
	$log->info("Running pipeline step #6: smrt bbinfer");
	
	$success = $self->_run_subcommand( ("bbinfer",  "-w", $workdir, "-s", $supermatrixfile, "-t", $classtreefile, "-o", $backbonefile, "-i", "examl", "-l", "bbinfer.log") );	
	die ("Pipeline step 'bbinfer' failed. Try to run 'smrt bbinfer' manually with '-v' option to find out what went wrong.") unless ( $success and  -e $backbonefile and -s $backbonefile);
	$log->info("Step #6 smrt bbinfer succeed");

	# Step 7: BBreroot
	$log->info("Running pipeline step #7: smrt bbreroot");
	if ( $outgroup ) {
		$success = $self->_run_subcommand( ("bbreroot",  "-w", $workdir, "-b", $backbonefile, "-t", $taxafile, "-o", $rerootedfile, "-g", $outgroup, "-l", "bbreroot.log") );	
	}
	else {
		$success = $self->_run_subcommand( ("bbreroot",  "-w", $workdir, "-b", $backbonefile, "-t", $taxafile, "-o", $rerootedfile, "-l", "bbreroot.log") );			
	}
	die ("Pipeline step 'bbreroot' failed. Try to run 'smrt bbreroot' manually with '-v' option to find out what went wrong.") unless ( $success and  -e $rerootedfile and -s $rerootedfile);
	$log->info("Step #7 smrt bbreroot succeed");
	
	# Step 8: BBcalibrate
	if ( my $fossilfile = $opt->fossiltable ) {
		$log->info("Running pipeline step #8: smrt bbcalibrate");
		$success = $self->_run_subcommand( ("bbcalibrate",  "-w", $workdir, "-s", $supermatrixfile, "-t", $rerootedfile, "-o", $chronogramfile, "-f", $fossilfile, "-l", "bbcalibrate.log") );	
		die ("Pipeline step 'bbcalibrate' failed. Try to run 'smrt bbcalibrate' manually with '-v' option to find out what went wrong.") unless ( $success and  -e $chronogramfile and -s $chronogramfile);
		$log->info("Step #8 smrt bbcalibrate succeed");					
	}
	else {
		$log->info("No fossil table given as argument. Skipping pipeline step #7: smrt bbcalibrate");
		$chronogramfile = $backbonefile;
	}
	
	# Step 9: BBdecompose
	$log->info("Running pipeline step #9: smrt bbdecompose");
	$success = $self->_run_subcommand( ("bbdecompose",  "-w", $workdir, "-b", $chronogramfile, "-c", $classtreefile, "-a", $alnfile, "-t", $taxafile, "-l", "bbdecompose.log") );	
	die ("Pipeline step 'bbdecompose' failed. Try to run 'smrt bbdecompose' manually with '-v' option to find out what went wrong.") unless ( $success );
	$log->info("Step #9 smrt bbdecompose succeed");
	
	
	# Step 8: Clademerge
	$log->info("Running pipeline step #10: smrt clademerge");
	$success = $self->_run_subcommand( ("clademerge",  "-w", $workdir, "-l", "clademerge.log") );	
	die ("Pipeline step 'clademerge' failed. Try to run 'smrt clademerge' manually with '-v' option to find out what went wrong.") unless $success;
	$log->info("Step #10 smrt clademerge succeed");
		
	# Step 9: Cladeinfer
	$log->info("Running pipeline step #11: smrt cladeinfer");
	$success = $self->_run_subcommand( ("cladeinfer",  "-w", $workdir, "-n", $beast_gens, "-l", "cladeinfer.log") );	
	die ("Pipeline step 'bbinfer' failed. Try to run 'smrt bbinfer' manually with '-v' option to find out what went wrong.") unless $success;
	$log->info("Step #11 smrt cladeinfer succeed");

	# Step 10: Cladegraft
	$log->info("Running pipeline step #12: smrt cladegraft");
	$success = $self->_run_subcommand( ("cladegraft",  "-w", $workdir, "-b", $chronogramfile, "-o", $finaltreefile, "-l", "cladegraft.log") );	
	die ("Pipeline step 'cladegraft' failed. Try to run 'smrt cladegraft' manually with '-v' option to find out what went wrong.") unless ( $success and  -e $finaltreefile and -s $finaltreefile);
	$log->info("Final step #12 smrt cladegraft succeed. Final file written to $finaltreefile");

	$log->info("DONE");
	
	return 1;
}

sub _run_subcommand {
	my ($self, @opt) = @_;
	local @ARGV = @opt;
	return Bio::SUPERSMART::App::smrt->run;	
}

1;