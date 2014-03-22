package Bio::Phylo::PhyLoTA::Service::TreeService;
use strict;
use warnings;
use Cwd;
use File::Temp 'tempfile';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service;
use base 'Bio::Phylo::PhyLoTA::Service';

my $config = Bio::Phylo::PhyLoTA::Config->new;
my $log = Bio::Phylo::PhyLoTA::Service->logger;	

# writes a tree description that spans provided taxa/clade to file
sub write_tree {

}

# it appears (obviously) that midpoint rooting is going to fail in some
# cases. How to automate this? One possibility is that the common tree might
# be of use.
sub reroot_tree {
	
}

# -infile, (optional: -burnin)
sub consense_trees {
	my ( $self, %args ) = @_;
	my $burnin = $args{'-burnin'} || $config->BURNIN;
	my $infile = $args{'-infile'} || die "Need -infile argument!";
	
	# we first need to count the absolute number of trees in the file
	# and multiply that by burnin to get the absolute number of trees to skip
	my $counter = 0;
	open my $fh, '<', $infile or die $!;
	while(<$fh>) {
		$counter++ if /^\s*tree\s/i;
	}	
	my $babs = int( $burnin * $counter );
	
	# create temporary outfile name
	my ( $outfh, $outfile ) = tempfile();
	close $outfh;
	
	# execute command
	my $tmpl = '%s -burnin %i %s %s 2> /dev/null';
	my $command = sprintf $tmpl, $config->TREEANNOTATOR_BIN, $babs, $infile, $outfile;
	system($command) and die "Error building consensus: $?";
	# return the resulting tree
	return parse_tree(
	        '-format' => 'nexus',
		'-file'   => $outfile,
		'-as_project' => 1,
	);
}

# grafts a clade tree into a backbone tree, returns backbone
sub graft_tree {
	my ( $self, $backbone, $clade ) = @_;
	my @names = map { $_->get_name } @{ $clade->get_terminals };
	$log->debug("@names");
	
	# retrieve the tips from the clade tree that also occur in the backbone, i.e. the 
	# exemplars, and locate their MRCA on the backbone
	my @exemplars;
	for my $name ( @names ) {
		if ( my $e = $backbone->get_by_name($name) ) {
			$log->info("found $name in backbone: ".$e->get_name);
			push @exemplars, $e;
		}
		else {
			$log->debug("$name is not in the backbone tree");
		}
	}
	$log->info("found ".scalar(@exemplars)." exemplars in backbone");
	my @copy = @exemplars; # XXX ???	
	my $bmrca = $backbone->get_mrca(\@copy);
	
	# find the exemplars in the clade tree and find *their* MRCA
	my @clade_exemplars = map { $clade->get_by_name($_->get_name) } @exemplars;
	my @ccopy = @clade_exemplars;
	my $cmrca = $clade->get_mrca(\@ccopy);
	$log->info("found equivalent ".scalar(@clade_exemplars)." exemplars in clade");
	
	# calculate the respective depths, scale the clade by the ratios of depths
	my $cmrca_depth = $cmrca->calc_max_path_to_tips;
	my $bmrca_depth = $bmrca->calc_max_path_to_tips;
	$clade->visit(sub{
		my $node = shift;
		my $bl = $node->get_branch_length || 0; # zero for root
		$node->set_branch_length( $bl * ( $bmrca_depth / $cmrca_depth ) );
		$log->debug("adjusting branch length for ".$node->get_internal_name);
	});
	$log->info("re-scaled clade by $bmrca_depth / $cmrca_depth");
	
	# calculate the depth of the root in the clade, adjust backbone depth accordingly
	my $crd  = $clade->get_root->calc_max_path_to_tips;
	my $diff = $bmrca_depth - $crd;
	$bmrca->set_branch_length( $bmrca->get_branch_length + $diff );
	$log->info("adjusted backbone MRCA depth by $bmrca_depth - $crd");
	
	# now graft!
	$bmrca->clear();
	$clade->visit(
		sub {
			my $node = shift;
			my $name = $node->get_internal_name;
			if ( my $p = $node->get_parent ) {
				if ( $p->is_root ) {
					$log->info("grafted $name onto backbone MRCA");
					$node->set_parent($bmrca);
				}
			}
			if ( $node->is_root ) {
				$log->info("replacing root $name with backbone MRCA");
			}
			else {
				$log->debug("inserting $name into backbone");
				$backbone->insert($node);
			}
		}
	);
	return $backbone;
}

# builds a tree using the configured tree inference method
sub infer_backbone_tree {
        
}

sub infer_clade_tree {
        
}

sub make_phylip_binary {
	my ( $self, $phylip, $binfilename, $parser, $work_dir) = @_;
	$log->info("going to make binary representation of $phylip => $binfilename");	
        $log->info("using parser $parser");	
	my ( $phylipvolume, $phylipdirectories, $phylipbase ) = File::Spec->splitpath( $phylip );
        my $curdir = getcwd;
	chdir $work_dir;        
	my @command = ( $parser, 
		'-m' => 'DNA', 
		'-s' => $phylipbase, 
		'-n' => $binfilename,
		'>'  => File::Spec->devnull,		
		'2>' => File::Spec->devnull,
	);
	my $string = join ' ', @command;
	$log->info("going to run '$string' inside " . $work_dir );
	my $a = system($string) and $log->warn("Couldn't create $binfilename: $?");        
	chdir $curdir;
	return "${binfilename}.binary";
}

sub make_phylip_from_matrix {
	my ( $self, $taxa, $phylipfile, @matrix ) = @_;
	
	# create phylip file for parser
	open my $phylipfh, '>', $phylipfile or die $!;
	my %nchar_for_matrix;
	my $ntax  = $taxa->get_ntax;
	my $nchar = 0;
	for my $m ( @matrix ) {
		my $mid = $m->get_id;
		$nchar_for_matrix{$mid} = $m->get_nchar;
		$nchar += $nchar_for_matrix{$mid};
	}
	print $phylipfh $ntax, ' ', $nchar, "\n";
	$taxa->visit(sub{
		my $t = shift;
		my @d = @{ $t->get_data };
		print $phylipfh $t->get_name, ' ';
		for my $m ( @matrix ) {
			my $mid = $m->get_id;
			my ($row) = grep { $_->get_matrix->get_id == $mid } @d;
			if ( $row ) {
				print $phylipfh $row->get_char;
			}
			else {
				print $phylipfh '?' x $nchar_for_matrix{$mid};
			}
		}
		print $phylipfh "\n";
	});	
	return $phylipfile;
}



# reads the supermatrix (phylip format) and returns the tip names from it
sub read_tipnames {
	my ($self, $supermatrix) = @_;
        my $logger = $self->logger;	
	my $ntax = 0;
	my $line = 0;
	my @result;	
        open my $fh, '<', $supermatrix or die $!;
	LINE: while(<$fh>) {
		chomp;
		my $word;                     
		if ( /^(\S+)/ ) {                        
			$word = $1;
                }
		if ( not $ntax ) {
			$ntax = $word;                        
			$logger->debug("$supermatrix has $ntax taxa");
                        next LINE;
                }
		if ( $word ) {
                        push @result, $word;
                };
		$logger->debug("adding taxon $word");
		last LINE if ++$line == $ntax;
	}
	return @result;
}


1;
