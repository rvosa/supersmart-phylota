#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::PhyLoTA::Service;

# process command line arguments
my ( $alignments, $product, $gene, $label );
my $verbosity = WARN;
GetOptions(
	'alignments=s' => \$alignments,
	'verbose+'     => \$verbosity,
	'product=s'    => \$product,
	'gene=s'       => \$gene,
	'label=s'      => \$label,
);

# instantiate helper objects
my $srv = Bio::Phylo::PhyLoTA::Service->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-class' => 'main',
	'-level' => $verbosity,
);

# read list of alignment files
my @alignments;
{
	$log->info("going to read list of alignments from $alignments");
	open my $fh, '<', $alignments or die $!;
	while(<$fh>) {
		chomp;
		push @alignments, $_;
	}
	close $fh;
}


# compose query
my %query;
if ( $product ) {
	$query{product} = $product;
}
if ( $gene ) {
	$query{gene} = $gene;
}

# iterate over alignments
for my $aln ( @alignments ) {
	$log->info("going to evaluate $aln");
	if ( $aln =~ /(\d+)\.\.fa$/ ) {
		my $gi = $1;
		$query{gi} = $gi;
		$log->info("seed GI is $gi");

		# this queries whether we even have an annotation for this seq
		$log->debug(Dumper(\%query));
		if ( my $feature = $srv->single_feature(\%query) ) {
			$log->info("found single feature that matches query");

			# need to have start and length, otherwise can't continue
			my $start  = $feature->codon_start;
			my $length = $feature->length;
			if ( defined $start and defined $length ) {
				$log->info("feature starts at $start and is $length bp");

				# now read the file
				my %fasta = simple_fasta(slurp_file($aln));
				$log->info("read file contents");

				# converts coordinates relative to unaligned sequence
				my @coords;
				my ($seq) = map { $fasta{$_} } grep { /gi\|$gi\|/ } keys %fasta;
				my @cols = split //, $seq;
				my $pos = 0;
				for my $i ( 0 .. $#cols ) {
					if ( $cols[$i] ne '-' and $cols[$i] ne '?' ) {
						$pos++;
					}
                                        if ( $pos >= $start and $pos <= ( $start + $length ) ) {
                                        	push @coords, $i;
                                        }
				}
				$log->info("computed coordinates @coords");

				# generate the new file name
				my $newfile = $aln;
				$newfile =~ s/\.\.fa$/.$label.fa/;
				$log->info("will write results to file $newfile");

				# now write the spliced contents
				open my $fh, '>', $newfile or die $!;
				for my $defline ( keys %fasta ) {
					print $fh ">$defline\n";
					my @seq = split //, $fasta{$defline};
					print $fh @seq[@coords], "\n";
				}
				close $fh;
				
			}
			else {
				$log->warn("no start or length for feature ".$feature->feature_id);
			}
		}
		else {
			$log->warn("no exact feature hit for sequence $gi");
			my @features = $srv->search_feature({ gi => $gi })->all;
			$log->info("sequence $gi has ".scalar(@features)." featuresÂ");
		}
	}
	else {
		$log->info("skipping file $aln");
	}
}


# reads entire file contents into a string
sub slurp_file {
	my $file = shift;
	open my $fh, '<', $file or die $!;
	my $result = do { local $/; <$fh> };
	close $fh;
	return $result;
}



# reads a FASTA string, returns a hash keyed on definition line (sans '>'
# prefix), value is concatenated seq
sub simple_fasta {
    my $string = shift;
    my @lines = split /\n/, $string;
    my %fasta;
    my $current;
    for my $line ( @lines ) {
        chomp $line;
        if ( $line =~ /^>(.+)/ ) {
            $current = $1;
            if ( exists $fasta{$current} ) {
                $fasta{$current} = '';
            }            
        }
        else {
            $fasta{$current} .= $line;
        }
    }
    return %fasta;
}
