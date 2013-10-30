#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use File::Basename;
use File::Path qw'make_path rmtree';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Domain::BigTree;
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;

# process command line arguments
my ( $dir, $replicates, $verbosity, $seqfile, $treefile, $outfile, $consense, $execute, $clean ) = ( '.', 100, WARN );
GetOptions(
	'workdir=s'    => \$dir,
	'replicates=i' => \$replicates,
	'seqfile=s'    => \$seqfile,
	'treefile=s'   => \$treefile,
	'outfile=s'    => \$outfile,
	'verbose+'     => \$verbosity,
	'consense'     => \$consense,
	'clean'        => \$clean,
	'execute'      => \$execute,
);

# instantiate helper objects
my $conf = Bio::Phylo::PhyLoTA::Config->new;
my $btr  = Bio::Phylo::PhyLoTA::Domain::BigTree->new;
my $mts  = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;
my $log  = Bio::Phylo::Util::Logger->new(
	'-class' => 'main',
	'-level' => $verbosity,
);

# parse fasta file
my %fasta = $mts->parse_fasta_file($seqfile);
$log->info("read FASTA data from $seqfile");

# compute number of characters
my ($nchar) = map { length($_) } values %fasta;
my %taxa = map { $_ => 1 } $mts->get_taxa_from_fasta(%fasta);
my $ntax = scalar keys %taxa;

# create file stem
my $base = basename($seqfile);
my $stem = $dir . '/' . $base;
make_path($dir) if not -d $dir;

#Êone-based
my @treefiles;
for my $i ( 1 .. $replicates ) {
	my %seen;
	my $bootstrap_file = $stem . '.rep' . $i . '.phy';
	
	# create bootstrap replicate if not exists
	if ( not -e $bootstrap_file ) {
		$log->info("going to write bootstrap replicate $i ($bootstrap_file)");
		
		# create an array of indices sampled with replacement
		my @indices;
		while (scalar(@indices) < $nchar ) {
			push @indices, int rand $nchar;
		}
		my @sorted = sort { $a <=> $b } @indices;
		
		# write the bootstrapped data set as PHYLIP
		open my $fh, '>', $bootstrap_file;
		print $fh "$ntax $nchar\n";
		for my $row ( keys %fasta ) {
			if ( $row =~ /taxon\|(\d+)/ ) {
				my $taxon = $1;
				if ( not $seen{$taxon}++ ) {
					my @seq = split //, $fasta{$row};
					my $seq = join '', @seq[@sorted];				
					print $fh $taxon, ' ' x ( 10 - length($taxon) ), $seq, "\n";
				}
				else {
					$log->warn("already seen taxon $taxon");
				}
			}
		}
		close $fh;
	}
	
	# run the search if requested
	if ( $execute ) {		
		$log->info("going to run tree search on $bootstrap_file");
		push @treefiles, $btr->build_tree(
			'-seq_file'  => $bootstrap_file,
			'-tree_file' => $treefile,
			'-work_dir'  => $dir,
		);
	}
}

# build consensus tree if requested
if ( $consense ) {
	$log->info("going to build consensus tree");
	
	# concatenate all tree files to an intree for consense
	{
		open my $fh, '>', "$dir/intree" or die $!;
		for my $file ( @treefiles ) {
			open my $treefile, '<', $file or die "Couldn't open $file: $!";
			my $newick = do { local $/; <$treefile> };
			print $fh $newick;
		}
		close $fh;
	}
	
	# run consense, have to change CWD for it and fool it into thinking
	# that a command-line user is interacting with it
	{
		my $cwd = getcwd;
		chdir $dir;
		#my ( $chld_out, $chld_in );
		#my $pid = open2( $chld_out, $chld_in, $conf->CONSENSE_BIN );
		#print $chld_in "r\n"; # tell consense that the trees are rooted
		#print $chld_in "y\n"; # tell consense to accept all other settings
		my $command = $conf->CONSENSE_BIN . ' >/dev/null';
		open my $consense, '|-', $command or die $!;
		print $consense "r\n";
		print $consense "y\n";
		chdir $cwd;
	}
	
	# now we have to re-scale the "branch lengths", i.e. bipartition counts
	my $rescaled;
	{
		open my $fh, '<', "$dir/outtree" or die $!;
		my $newick = '';
		while(<$fh>) {
			chomp;
			$newick .= $_;
		}
		my @parts = split /:/, $newick;
		for my $i ( 0 .. $#parts ) {
			if ( $parts[$i] =~ /^(\d+\.?\d*)/ ) {
				my $length = $1;
				my $replace = $length / $replicates;
				$parts[$i] =~ s/^$length/$replace/;
			}
		}
		$rescaled = join ':', @parts;
		$rescaled =~ s/(\d+):\d+\.?\d*/$1/g; # strip the values at the tips
	}
	
	# write to the outfile
	{
		open my $fh, '>', $outfile or die $!;
		print $fh $rescaled;
		$log->info("tree written to $outfile");
	}
}

if ( $clean ) {
	$log->info("cleaning up working directory $dir");
	rmtree($dir);
}


