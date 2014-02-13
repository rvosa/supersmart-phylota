#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';

=head1 NAME

join_alignments.pl - joins alignments by taxon id

=head1 SYNOPSYS

 $ join_alignments.pl --list=<file> --format=<nexus|fasta|phylip> > <outfile>

=head1 DESCRIPTION

Given a text file that contains a list of file locations (one per line), produces an
output that joins these files on their taxon identifiers. The output is written to
STDOUT in either FASTA, NEXUS or PHYLIP format.

=cut

# process command line arguments
my $verbosity = WARN;
my ( $list, $format );
GetOptions(
	'list=s'   => \$list,
	'verbose+' => \$verbosity,
	'format=s' => \$format,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-class' => 'main',
	'-level' => $verbosity,
);

# read list of files
my @list;
{
	my $fh; # file handle
	
	# may also read from STDIN, this so that we can pipe 
	if ( $list eq '-' ) {
		$fh = \*STDIN;
		$log->debug("going to read file names from STDIN");
	}
	else {
		open $fh, '<', $list or die $!;
		$log->debug("going to read file names from $list");
	}
	
	# read lines into array
	@list = <$fh>;
	chomp @list;
}

# populate supermatrix
my %supermatrix;
my $nchar = 0;
my %charset;
for my $file ( @list ) {

	# read alignment file, return hash where key is NCBI taxon ID,
	# value is aligned sequence
	my %matrix = parse_matrix($file);
	my $matrix_nchar = length $matrix{ ( keys %matrix )[0] };
	$charset{$file} = [ $nchar + 1 => $nchar + $matrix_nchar ];
	
	# concat current data to supermatrix
	for my $row ( keys %matrix ) {	
	
		# if current data has row not yet seen, pad with '?'	
		if ( not exists $supermatrix{$row} ) {
			$supermatrix{$row} = '?' x $nchar;
		}		
		$supermatrix{$row} .= $matrix{$row};
	}
	
	# add missing for supermatrix rows not in matrix
	for my $row ( keys %supermatrix ) {
		if ( not exists $matrix{$row} ) {
			$supermatrix{$row} .= '?' x $matrix_nchar;
		}
	}
	
	# add current matrix's length to overall length
	$nchar += $matrix_nchar;
}

# write results
if ( $format =~ /nexus/i ) {
	write_nexus( $nchar, \%supermatrix, \%charset );
}
elsif ( $format =~ /fasta/i ) {
	for my $key ( keys %supermatrix ) {
		print '>taxon|', $key, "\n", $supermatrix{$key}, "\n";
	}
}
else {
	write_phylip( $nchar, \%supermatrix );
}

# writes a phylip representation of the supermatrix
sub write_phylip {
	my ( $nchar, $matrix ) = @_;
	my @taxa = sort { $a <=> $b } keys %{ $matrix };
	print scalar(@taxa), ' ', $nchar, "\n";
	print $_, ' ', $matrix->{$_}, "\n" for @taxa;
}

# writes a nexus representation of the supermatrix, with character sets
sub write_nexus {
	my ( $nchar, $matrix, $charset ) = @_;
	my @taxa = sort { $a <=> $b } keys %{ $matrix };
	
	# write taxa block
	print "#NEXUS\nBEGIN TAXA;\n";
	print "\tDIMENSIONS NTAX=".scalar(@taxa).";\n";
	print "\tTAXLABELS\n";
	print "\t\tt$_\n" for @taxa;
	print "\t;\nEND;\n";
	
	# write characters block
	print "BEGIN CHARACTERS;\n";
	print "\tDIMENSIONS NCHAR=${nchar};\n\tFORMAT DATATYPE=DNA GAP=- MISSING=?;\n";
	print "\tMATRIX\n";
	print "\t\tt$_\t", $matrix->{$_}, "\n" for @taxa;
	print "\t;\nEND;\n";
	
	# write sets block
	print "BEGIN SETS;\n";
	for my $set ( sort { $charset->{$a}->[0] <=> $charset->{$b}->[0] } keys %{ $charset } ) {
		my $start = $charset->{$set}->[0];
		my $end   = $charset->{$set}->[1];
		print "\tCHARSET '$set' = ${start}-${end};\n";
	}
	print "END;\n";
}

sub parse_matrix {
	my $file = shift;
	my ( %matrix, $current );
	
	# open file handle 
	open my $fh, '<', $file or die $!;
	$log->info("going to read sequences from FASTA file $file");
	
	# read over the file handle
	while(<$fh>) {
		chomp; # strip line ending
		
		# this matches the FASTA definition line, we capture the taxon ID
		if ( />.*taxon\|(\d+)/ ) {
			$current = $1;
			$log->debug("found taxon ID $current");
			
			# we have already seen this ID, now we're seeing it for the
			# second (or more) time.
			if ( $matrix{$current} ) {
				$log->debug("already seen $current, starting new empty string");
				push @{ $matrix{$current} }, '';
			}
			
			# this is the first time we see the taxon ID
			else {
				$log->debug("not yet seet $current, initializing array with an empty string");
				$matrix{$current} = [ '' ];
			}
		}
		else {			
			s/\s//g;
			$matrix{$current}->[-1] .= $_;
		}
	}
	
	# now we're going to reduce the sequence set to use only the longest
	# sequence in the set for that taxon ID
	$log->info("going to select longest sequence from among within-species sequences");
	for my $taxon_id ( keys %matrix ) {
		my @sequences = @{ $matrix{$taxon_id} };
		
		# we sort in ascending order from fewer to more gaps
		my @sorted = sort { $a =~ tr/-/-/ <=> $b =~ tr/-/-/ } @sequences;
		$matrix{$taxon_id} = $sorted[0];
	}
	return %matrix;
}
