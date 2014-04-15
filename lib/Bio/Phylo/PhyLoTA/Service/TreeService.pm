package Bio::Phylo::PhyLoTA::Service::TreeService;
use strict;
use warnings;
use Cwd;
use File::Temp 'tempfile';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service;
use base 'Bio::Phylo::PhyLoTA::Service';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use Storable 'dclone';

my $config = Bio::Phylo::PhyLoTA::Config->new;
my $log = Bio::Phylo::PhyLoTA::Service->logger;	

# writes a tree description that spans provided taxa/clade to file
sub write_tree {

}

# re-roots a given backbone tree such that the number of monophyletic genera is maximized 
sub reroot_tree {       
	my ($self, $tree, @records) = @_;
        if (! @records) {
                die "Need both, tree and species table!";
        }
        $log->info("rerooting tree");
        my $nw = $tree->to_newick;

        # iterate over all nodes, re-root the tree at the branch below each
        # node and check how many clades we will get. Take the tree that will 
        # have the most clades.
        my $max_clades = scalar ( $self->_make_clade_species_sets( $tree, @records ) );
        my $num_internals = scalar @{ $tree->get_internals };
        my $result = $tree;

        foreach my $i ( 0..$num_internals-1 ) {
                # Get a fresh curr_tree, because node order could be messed up
                # by re-rooting and not all internal branches would be considered.
                # Reading from newick is a workaround since $tree->clone somehow takes forever
                my $curr_tree = parse_tree(
                        '-format'     => 'newick',
                        '-string'     => $nw,
                    );                                
                my $node = @{ $curr_tree->get_internals }[$i];                
                $curr_tree->reroot($node);
                my @set = $self->_make_clade_species_sets( $curr_tree, @records );
                if ( scalar (@set) > $max_clades ){
                        $max_clades = scalar (@set);
                        $result = $curr_tree;
                }
        }
        $log->info("rerooted tree has $max_clades monophyletic clades");
        return $result;
}

# Changes NCBI taxon identifiers in a given tree to taxon names
sub remap {
        my ($self, $tree) = @_;
        $tree->visit(sub{
                my $n = shift;
                if ( $n->is_terminal ) {
                        my $id = $n->get_name;
                        $n->set_name( $self->find_node($id)->taxon_name );
                }
                     });
        return $tree;
        
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
        my $str = $clade->to_newick;
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
			$log->debug("$supermatrix has $ntax taxa");
                        next LINE;
                }
		if ( $word ) {
                        push @result, $word;
                };
		$log->debug("adding taxon $word");
		last LINE if ++$line == $ntax;
	}
	return @result;
}



## as
sub _make_clade_species_sets {
        my ($self, $tree, @records) = @_;
        if (! @records) {
                die "Need both, tree and species table!";
        }                
        my $mt = 'Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa';

        # traverse tree, retrieve monophyletic genera, which we will compare to all
        $log->info("going to identify putatively monophyletic genera");
        my ( %monophyletic, %all );
        $tree->visit_depth_first(
                '-post' => sub {
                        my $node = shift;
                        
                        # for the tips, we keep a running tally of all seen genera. if a genus 
                        # occurs on more than one tip we need to know if the genus is monophyletic
                        if ( $node->is_terminal ) {
                                my $id = $node->get_name;
                                my ($genus) = map { $_->{'genus'} } grep { $_->{'species'} == $id } @records;
                                $node->set_generic( 'species' => [ $id ], 'genera' => { $genus => 1 } );
                                
                                # in the end, this will be either 1 (for monotypic genera), or 2 (for
                                # exemplar genera).
                                $all{$genus}++;
                        }
                        else {
                                
                                # here we simply carry over the species and genera from the
                                # children to the focal node. in the second pass we will use
                                # this to find the shallowest node that subtends all species
                                # of a paraphyletic genus
                                my ( @s, %g );						
                                for my $child ( @{ $node->get_children } ) {
                                        push @s, @{ $child->get_generic('species') };
                                        my %genus = %{ $child->get_generic('genera') };
                                        $g{$_} += $genus{$_} for keys %genus;
                                }
                                $node->set_generic( 'species' => \@s, 'genera' => \%g );
                                
                                # the node subtends two species that all belong to the same genus. 
                                # hence, the node is monophyletic.
                                if ( scalar(@s) >= 2 and scalar keys %g == 1 ) {
                                        $monophyletic{$_} += $g{$_} for keys %g;
                                }
                        }
                }
            );        
        
        # all the genera that aren't putatively monophyletic but that do have two
        # members are therefore paraphyletic. in the second pass we lump all the
        # putatively monophyletic genera that nest inside paraphyletic ones within
        # the mixed set of mono/para.
        $log->info("going to reconstruct paraphyly");
        my %paraphyletic = map { $_ => 1 } grep { !$monophyletic{$_} && $all{$_} == 2 } keys %all;
        my @genera;
        $tree->visit_depth_first(
                '-post' => sub {
                        my $node = shift;
                        if ( $node->is_internal ) {
                                my %g = %{ $node->get_generic('genera') };
                                
                                # the node is a paraphyletic mrca if it is the shallowest
                                # node where a paraphyletic genus occurs and where it
                                # subtends the two exemplars from that genus
                                my $is_para;
                                for my $genus ( keys %g ) {
                                        if ( $paraphyletic{$genus} and $g{$genus} == 2 ) {
                                                $is_para++;
                                        }
                                }
                                
                                # if the node is paraphyletic we store ALL its subtended
                                # genera, removing them from the set of putative monophyletic
                                # genera as well as from the paraphyletic ones. we remove
                                # from the former set so that after this traversal we don't
                                # put monophyletic genera that nest inside paraphyletic ones
                                # in a separate set, and we remove the latter set so that
                                # deeper nodes that also subtend the paraphyletic genera
                                # don't trigger processing.
                                if ( $is_para ) {
                                        my @g = keys %g;
                                        push @genera, \@g;
                                        delete @paraphyletic{@g};
                                        delete @monophyletic{@g};
                                }	
                        }	
                }
            );
        # all remaining ones become their own set
        push @genera, [ $_ ] for keys %monophyletic;
        # now resolve the nesting of paraphyletic genera
        my %index;
        for my $i ( 0 .. $#genera ) {
                $index{$_} = $i for @{ $genera[$i] };                
        }
        
        # make species sets
        my @set;
        for my $i ( keys %{ { map { $_ => 1 } values %index } } ) {
                my @g = @{ $genera[$i] };
                my @s = map { $mt->get_species_for_taxon( 'genus' => $_, @records ) } @g;
                push @set, \@s if scalar(@s) > 2;
        }
        
        return @set;        
}


1;
