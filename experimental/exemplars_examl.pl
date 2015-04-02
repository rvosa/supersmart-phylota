

sub get_exemplars_for_genus {
    my $genus = shift;
    my @records = @_;
    print "Genus : $genus\n";
    use Bio::Phylo::Util::Logger ':levels';
    my $logger =  Bio::Phylo::Util::Logger->new(
	'-level' => DEBUG,
	'-class' => [':Bio::SUPERSMART::App::smrt::Command::BBmerge']);
    my $config = Bio::Phylo::PhyLoTA::Config->new;
    # get species for that genus
    my $mt = Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa->new;    
    my @species = $mt->query_taxa_table([$genus], ['species'], @records);
    system("mkdir $genus");
    
    # now read the list of alignments
    my $alnfile = './aligned-local.txt';
    my @alignments;
    $logger->info("going to read list of alignments $alnfile");
    open my $fh, '<', $alnfile or die $!;
    while(<$fh>) {
	chomp;
	push @alignments, $_ if /\S/ and -e $_;
    }
    # iterate over alignments, extract the sequences for our species and   
    # write the new alignment into the genus directory
    for my $aln (@alignments) {
	$logger->info("Processing alignment $aln");
	my %fasta = $mt->parse_fasta_file($aln);
	open my $fh, '<', $aln or die $!;
	my $fastastr = join('', <$fh>);
	close $fh;
	my $dist  = $mt->calc_mean_distance($fastastr);
	my $nchar = $mt->get_nchar(%fasta);
	my $mdens = $config->CLADE_MIN_DENSITY;
	my %seq;
	my $distinct = 0;
	if ( $dist <= $config->CLADE_MAX_DISTANCE ) {
	    for my $s ( @species ) {
		my @s = grep { /taxon\|$s[^\d]/ } keys %fasta;
		if ( @s ) {
		    $seq{$s} = \@s;
		    $distinct++;
		}
	    }
	    if ( ($distinct/scalar(@species)) >= $mdens && $distinct > 2 ) {
		$logger->info("$aln is informative and dense enough");
		open my $fh, '>', "./$genus/$aln" or die $!;
		for my $s ( @species ) {
		    if ( $seq{$s} ) {
			for my $j ( 0 .. $#{ $seq{$s} } ) {
			    my $seq = $seq{$s}->[$j];
			    print $fh '>', $seq, "\n", $fasta{$seq}, "\n";
			}
		    }
		}
		close $fh;
	    }
	    else {
		my $dens = sprintf ( "%.2f", $distinct/scalar(@species));
		$logger->info("$aln is not informative or dense enough: density " 
			      . $dens . " < $mdens, number of distinct taxa: $distinct");
	    }
	} 
	else {
	    $logger->info("Alignment $aln too divergent!!!");
	}
    }        
}
