package Bio::SUPERSMART::App::smrt::Command::Cladeinfer;

use strict;
use warnings;

use File::Copy;
use File::Temp 'tempfile';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service::ParallelService;
use Bio::Phylo::Util::Exceptions 'throw';
use Bio::Tools::Run::Phylo::StarBEAST;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: infers a species tree for each individual clade using *BEAST.

=head1 NAME

Cladeinfer.pm - For each decomposed clade, infers the species tree. This is 
done using the multispecies, multilocus coalescent as implemented in *BEAST 
on the basis of the merged set of alignments for that clade.

=head1 SYNOPSYS

smrt cladeinfer [-h ] [-v ] [-w <dir>] [-n <int>] [-s <int>] [-l <int>] [-f ] \
  [--rebuild] [--append --burnin=0.1]

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

=back

=cut


sub options {
    my ($self, $opt, $args) = @_;
    my $config = Bio::Phylo::PhyLoTA::Config->new;
    return (
        [ "ngens|n=i", "number of generations in *BEAST, defaults to 100_000 (strictly for testing!)", { arg => "int", default => 100_000 } ],
        [ "sfreq|s=i", "sampling frequency, defaults to 1000", { arg => "int", default => 1000 } ],
        [ "lfreq|l=i", "logging frequency, defaults to 1000", { arg => "int", default => 1000 } ],
        [ "file|f=s", "file (nexml format) to start a single inference from", { arg=>"file" } ],
        [ "outfile|o=s", "location of output directory", { arg => "location", default => "cladeinfer_out.txt"} ],
        [ "rebuild|x", "rebuild existing *BEAST XML files", {}],
        [ "append|a", "append trees to existing tree file", {}],
        [ "burnin|b=f", "burnin to omit when appending trees", { arg => "float", default => $config->BURNIN } ],
    );
}

sub validate {
    my ($self, $opt, $args) = @_;

    if ( my $file = $opt->file ) {
        $self->usage_error("file $file does not exist") unless -e $file;
        $self->usage_error("file $file is empty") unless -s $file;
    }
    if ( (my $ngens = $opt->ngens) <  30_000_000 ) {
        $self->logger->warn("chain length $ngens seems very low, are you just testing?");
    }
    if ( $opt->append ) {    	
    	if ( $opt->burnin < 0 or $opt->burnin > 1 ) {
    		$self->usage_error("burnin should be number between 0 and 1");
    	}
    	if ( not $opt->rebuild ) {
    		$self->usage_error("--append requires the --rebuild flag");
    	}
    }
}

=item make_suffix

Generates a file suffix to append to tree and log file

=cut

sub make_suffix {
	my ( $self, $stem, $append ) = @_;
	my $suffix = '';
	my $i = 1;
	if ( $append ) {
		$i++ while -e "${stem}.nex." . $i;
		$suffix = '.' . $i;
	}
	return $suffix;	
}

=item append_logs

Appends new log data (trees) to existing files

=cut

sub append_logs {
	my ( $self, %args ) = @_;
	
	# closure to perform burnin generically
	my $burner = sub {
		my ( $burnin, @list ) = @_;
		my $index = int( $burnin * scalar @list );
		return @list[$index..$#list];
	};
	
	# append tree files
	my ( $source, $target ) = @{ $args{'trees'} };
	if ( -e $source and -e $target and $source ne $target ) {	
		$self->append_trees( $source, $target, $burner, $args{'burnin'} );
	}
	
	# XXX log data wasn't quite working, maybe need to update column 1?
}

=item append_trees

Appends new trees to existing file

=cut

sub append_trees {
	my ( $self, $source, $target, $burner, $burnin ) = @_;
	
	# parse out just the trees
	my @trees;
	open my $sfh, '<', $source or die $!;
	while(<$sfh>) {
		push @trees, $_ if /^tree/;
	}
	
	# write the target file to temp
	my ( $fh, $filename ) = tempfile();		
	my $tree = undef;
	open my $tfh, '<', $target or die $!;
	while(<$tfh>) {
		$tree = 1 if /^tree/;
		if ( $tree ) {
			print $fh $_ if /^tree/;
		}
		else {
			print $fh $_;
		}		
	}
	close $tfh;
	
	# write the trees to temp and close the nexus
	print $fh $burner->( $burnin => @trees ), 'End;';
	close $fh;
	
	# clean up
	copy($filename,$target) or die "Append failed: $!";
	unlink $filename;	
	unlink $source;
}

sub run {
    my ($self, $opt, $args) = @_;

    # collect command-line arguments
    my $ngens = $opt->ngens;
    my $sfreq = $opt->sfreq;
    my $lfreq = $opt->lfreq;
    my $file  = $opt->file;
    my $outfile = $self->outfile;   
    my $workdir = $self->workdir;   
    my $rebuild = $opt->rebuild;    
        
    # instantiate helper objects
    my $config = Bio::Phylo::PhyLoTA::Config->new;
    my $beast  = Bio::Tools::Run::Phylo::StarBEAST->new;
    my $logger = $self->logger;

    # configure beast
    $logger->info("setting beast template file to " . $config->BEAST_TEMPLATE_FILE);
    $beast->template($config->BEAST_TEMPLATE_FILE);

    $logger->info( "setting beast executable to " . $config->BEAST_BIN );
    $beast->executable( $config->BEAST_BIN );
    
    $logger->info("setting chain length to $ngens");
    $beast->chain_length($ngens);
    
    $logger->info("setting sampling frequency to $sfreq");
    $beast->sample_freq($sfreq);
    
    $logger->info("setting logging frequency to $lfreq");
    $beast->log_freq($lfreq);

    $logger->info("setting seed from $config: ".$config->RANDOM_SEED );
    $beast->seed($config->RANDOM_SEED);
    
    # XXX these should be configurable from phylota.ini
    $beast->beagle_SSE(1);
    $beast->beagle_CPU(1);
    $beast->beagle_instances(1);
    
    # overwrite any previously existing output files, $rebuild applies
    # to any previously existing input BEAST XML files, which may 
    # have been edited
    $beast->overwrite(1);
    $beast->rebuild($rebuild);
    
    # run either one file or a directory
    if ( $file and -e $file ) {
        (my $stem = $file) =~ s/\..+//g; 
        
        # generate optional suffix, for re-runs
        my $suffix = $self->make_suffix( $stem, $opt->append );
        if ( $suffix ) {
        	$logger->info("Will append '$suffix' to output files");
        }
        
        # set output file names
        $logger->info("Starting inference for single clade tree");
        $beast->outfile_name( "${stem}.nex${suffix}" );
        $beast->logfile_name( "${stem}.log${suffix}" );
        
        # set input file
        $logger->info("Setting beast input file name to ${stem}-beast-in.xml");
        $beast->beastfile_name( "${stem}-beast-in.xml" );
        $beast->run( $file );
        $logger->info("done. trees are in ${file}.nex, BEAST log in ${file}.log");
        
		# concatenate 
		if ( $opt->append ) {
			$self->append_logs(
				'trees'  => [ "${stem}.nex${suffix}" => "${stem}.nex" ],
				'params' => [ "${stem}.log${suffix}" => "${stem}.log" ],
				'burnin' => $opt->burnin,
			);
		}        
    }
    
    else {
      
        # iterate over entries in work dir
        my @cladedirs;
        opendir my $dh, $workdir or die $!;
        while( my $entry = readdir $dh ) {
        
            # peruse directories named cladeXXX
            if ( $entry =~ /clade\d+/ && -d "${workdir}/${entry}" ) {
                push @cladedirs, $entry;
            }
        }        
        die "no clade directories found in workdir $workdir" if scalar(@cladedirs) == 0;

        # infer clades in parallel mode
        pmap {
            my ($clade) = @_;

            # this should be a nexml file with one taxa block and
            # multiple matrices
            my $stem = "${workdir}/${clade}/${clade}";
            my $file = "${stem}.xml";
            if ( -e $file ) {

                # nexml file exists
                my $suffix = $self->make_suffix( $stem, $opt->append );
                $beast->outfile_name( "${stem}.nex${suffix}" );
                $beast->logfile_name( "${stem}.log${suffix}" );
                $beast->beastfile_name( "${stem}-beast-in.xml" );
                $beast->run( $file );
                
                # concatenate 
                if ( $opt->append ) {
					$self->append_logs(
						'trees'  => [ "${stem}.nex${suffix}" => "${stem}.nex" ],
						'params' => [ "${stem}.log${suffix}" => "${stem}.log" ],
						'burnin' => $opt->burnin,
					);
                }
                
                my $tmpl = 'done with %s. trees are in %s.nex, log is in %s.log';
                $logger->info(sprintf $tmpl, $clade, $stem, $stem);             
            }
            else {
                $logger->warn("inconsistent directory structure, missing: $file");
            }
        } @cladedirs;
    }
    open my $fh, '>', $outfile;
    print $fh "Cladeinfer done\n";
    close $fh;
    
    $logger->info("DONE."); 
}

1;
