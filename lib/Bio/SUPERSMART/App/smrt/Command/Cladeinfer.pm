package Bio::SUPERSMART::App::smrt::Command::Cladeinfer;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Config;
use Bio::Tools::Run::Phylo::StarBEAST;

use Bio::Phylo::PhyLoTA::Service::ParallelService 'pthreads'; 

use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: infers a species tree for each individual clade using *BEAST.

=head1 NAME

Cladeinfer.pm - For each decomposed clade, infers the species tree. This is 
done using the multispecies, multilocus coalescent as implemented in *BEAST 
on the basis of the merged set of alignments for that clade.

=head1 SYNOPSYS

smrt cladeinfer [-h ] [-v ] [-w <dir>] [-n <int>] [-s <int>] [-l <int>] [-f ] 

=head1 DESCRIPTION

Infers a species-level phylogeny using the multi-species, multi-locus coalescent as
implemented in *BEAST. Given a directory, traverses it, looks for subdirectories
and files that match the pattern clade*.xml>. These must be NeXML files
with one or more alignments that link to the same taxa block in them. Given a file
argument, analyses a single file. *BEAST is run for each clade with one chain. Here we 
make use of parallel processing by simultaneously running *BEAST for different clades
on different processors. 

=over

=item ngens

Specifies the number of generations. The default, which is strictly for testing, is 100.
Instead this should be something that is more on the order of 30_000_000 if not much
more.

=item sfreq
    
Specifies the sampling frequency. The default, which is strictly for testing, is 100.
Instead this should be something that is more on the order of 300_000 if not much
more. For example, this could be 0.1% of ngens so that you end up with a thousand trees.

=item lfreq

Species the logging frequency. The default, which is strictly for testing, is 100.
Instead this should be something that is more on the order of 300_000 if not much
more. For example, this could be 0.1% of ngens so that you end up with a thousand logged
samples.

=cut


sub options {
	my ($self, $opt, $args) = @_;
	return (
		[ "ngens|n=i", "number of generations in *BEAST, defaults to 100_000 (strictly for testing!)", { arg=>"int", default=>100_000} ],
		[ "sfreq|s=i", "sampling frequency, defaults to ngens/100", { arg=>"int"} ],
		[ "lfreq|l=i", "logging frequency, defaults to ngens/100", { arg=>"int"} ],
		[ "file|f=s", "file (nexml format) to start a single inference from" ],
    );
}

sub validate {
	my ($self, $opt, $args) = @_;		

	if ( my $file = $opt->file ) {
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
	if ( (my $ngens = $opt->ngens) <  30_000_000 ) {
		$self->logger->warn("chain length $ngens seem very low, are you just testing?");		
	}
	if ( my $sfreq = $opt->sfreq <  300_000 ) {
		$self->logger->warn("sampling frequency $sfreq seem very low, are you just testing?");		
	}
	if ( my $lfreq = $opt->lfreq <  300_000 ) {
		$self->logger->warn("logging frequency $lfreq seems very low, are you just testing?");		
	}		
}

sub run {
	my ($self, $opt, $args) = @_;

   	# collect command-line arguments
    my $ngens = $opt->ngens;
    my $sfreq = $opt->sfreq ? $opt->sfreq : $ngens/100;
    my $lfreq = $opt->lfreq ? $opt->lfreq : $ngens/100;
	my  $file = $opt->file;
	(my $workdir = $opt->workdir) =~ s/\/$//g;
	
		
    # instantiate helper objects
    my $config = Bio::Phylo::PhyLoTA::Config->new;
    my $beast  = Bio::Tools::Run::Phylo::StarBEAST->new;
	my $logger = $self->logger;

	# overwrite existing output files
	my $overwrite = 1;

	# configure beast
	$logger->info( "setting beast executable to " . $config->BEAST_BIN );
	$beast->executable( $config->BEAST_BIN );
	
	$logger->info("setting chain length to $ngens");
	$beast->chain_length($ngens);
	
	$logger->info("setting sampling frequency to $sfreq");
	$beast->sample_freq($sfreq);
	if ( $sfreq < 300_000 ) {
		$logger->warn("sampling frequency $sfreq seem very low, are you just testing?");
	}
	
	$logger->info("setting logging frequency to $lfreq");
	$beast->log_freq($lfreq);
	if ( $lfreq < 300_000 ) {
		$logger->warn("logging frequency $lfreq seems very low, are you just testing?");
	}
	
	# XXX these should be configurable from phylota.ini
	$beast->beagle_SSE(1);
	$beast->beagle_CPU(1);
	$beast->beagle_instances(1);
	
	# overwrite any previously existing files
	$beast->overwrite($overwrite);
	
	# run either one file or a directory
	if ( $file and -e $file ) {
		$beast->outfile_name( "${file}.nex" );
		$beast->logfile_name( "${file}.log" );
		$logger->info("Setting beast input file name to ${file}.bxml");
		$beast->beastfile_name( "${file}-beast-in.xml" );
		$beast->run( $file );
		$logger->info("done. trees are in ${file}.nex, BEAST log in ${file}.log");
	}
	
	else {	
		# iterate over entries in work dir
		opendir my $dh, $workdir or die $!;
	    my @cladedirs;
		
		sequential { 
				while( my $entry = readdir $dh ) {
				# peruse directories named cladeXXX
					if ( $entry =~ /clade\d+/ && -d "${workdir}/${entry}" ) {
			                        push   @cladedirs, $entry;
			        }
		    }        
		};
		die("no clade directories found in workdir $workdir") if scalar(@cladedirs) == 0;
			
    # infer clades in parallel mode
    pmap {
                my $clade = $_;
                # this should be a nexml file with one taxa block and
                # multiple matrices                
                my $stem = "${workdir}/${clade}/${clade}";
                my $file = "${stem}.xml";
                if ( -e $file ) {
                        $beast->outfile_name( "${stem}.nex" );
                        $beast->logfile_name( "${stem}.log" );
						$beast->beastfile_name( "${stem}-beast-in.xml" );
                        $beast->run( $file );
                        my $tmpl = 'done with %s. trees are in %s.nex, log is in %s.log';
                        $logger->info(sprintf $tmpl, $clade, $stem, $stem);				
                }
                else {
                        $logger->warn("inconsistent directory structure, missing: $file");
                }
        } @cladedirs;
	}
	
	$logger->info("DONE.");	
}


1;