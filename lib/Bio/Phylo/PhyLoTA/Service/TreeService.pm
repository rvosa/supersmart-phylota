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
use List::MoreUtils 'uniq';
use List::Util 'min';

my $config = Bio::Phylo::PhyLoTA::Config->new;
my $log = Bio::Phylo::PhyLoTA::Service->logger->new;	

sub reroot_tree {
	my ( $self, $tree, $taxatable, $levels ) = @_;
	
	my @records = @{$taxatable};
	
	# taxonomic ranks to consider for determining paraphyly
	my @levels = @{$levels};
	
	# store the scores (number of paraphyletic species per rerooting) for each level and node index
	# at which we reroot
	my %scores = ();
	@scores{ @levels } = ();

	# store the tree objects for the rerooting at each node
	my %rerooted_trees = ();

	my $num_internals = scalar @{ $tree->get_internals };
	my $newi          = $tree->to_newick;

	# iterate over all possible rerootings
	for my $i (0 .. $num_internals - 1 ) {

		# get fresh tree, so the node order etc. won't be messed up
		my $current_tree = parse_tree(
			'-string' => $newi,
			'-format' => 'newick',
		);

		my @internals = @{ $current_tree->get_internals };

		# reroot the tree
		my $node = $internals[$i];
		$log->debug("rerooting tree at internal node # " . ($i+1) . "/" . $num_internals);

		$node->set_root_below;

		# store the tree for later
		$rerooted_trees{$i} = $current_tree;

		# calculate the 'scores' : count the amount of paraphyly at each taxonomic level
		foreach my $level (@levels) {
			my $para_count =
				$self->_count_paraphyletic_species( $current_tree, $level, @records );
			$log->debug("Paraphyletic species (counted from level $level) : $para_count");
			$scores{$level}{$i} = $para_count;
		}
	}

    # traverse all scores and just keep the node indices for the min scores for each level
	my %best_node_idx = ();
	foreach my $level (@levels) {
		$best_node_idx{$level} = ();
		my $minscore = min( values %{$scores{$level}} );
		foreach my $k ( keys %{$scores{$level}} ) {
			if ( $scores{$level}{$k} == $minscore ) {
				push @{ $best_node_idx{$level} }, $k;
			}
		}
	}

	# now we take the intersection of all indices for all the levels
	my @best_indices =
	  @{ $best_node_idx{ $levels[0] } };   #initialize with set from first level
		
	my $count = @levels;
	for my $ii ( 1 .. $count - 1 ) {
		@best_indices =
		  _intersect( \@best_indices, $best_node_idx{ $levels[$ii] }  );
	}

	if ( scalar @best_indices > 1 ) {
		$log->warn ("Found more than one optimal rerooted trees. Returning the first one");
	}
	if ( scalar @best_indices == 0 ) {
		$log->warn ("Found no optimal rerooted tree. Returning the original one");
		return $tree;
	}
	return $rerooted_trees{ $best_indices[0] };
}

sub _count_paraphyletic_species {
	my ( $self, $current_tree, $level, @records ) = @_;

	# get all distinct members of the NCBI tree at the given level
	my @members = uniq( map { $_->{$level} } @records );

	# sum up the amount of paraphyletic species for all members at that level
	my $para_count = 0;

	for my $m (@members) {
		
	   # get the subset of terminals in the NCBI tree that belong to this family
		my @ncbi_terminal_names =
		  map { $_->{'species'} } grep { $_->{$level} eq $m } @records;

		# get node objects for species in our tree
		my @nodes = map{ $current_tree->get_by_name( $_ )  } @ncbi_terminal_names;

		# get most recent common ancestor of the species in this family in our tree
		my $mrca = $current_tree->get_mrca( \@nodes );
		if ($mrca) {

			# get descendants of the mrca
			my @terminals = @{ $mrca->get_terminals };
			my @terminal_names = map { $_->get_name } @terminals;

			# count how many terminals of the mrca in our tree are in the
			#  terminals in the NCBI tree, paraphyletic species are then the
			#  ones that are descendents in our tree but not in the NCBI tree
			my @matching = _intersect( \@ncbi_terminal_names, \@terminal_names );
			my @paraspecies = _array_minus( \@terminal_names, \@ncbi_terminal_names );
			$para_count += scalar @paraspecies;
		}
	}
	return $para_count;
}

# get intersection between elements of two array references
sub _intersect {
	my ($a1, $a2) = @_;
	my @arr1 = @{$a1};
	my @arr2 = @{$a2};
	my %arr1 = map{$_=>1} @arr1;
	my %arr2 = map{$_=>1} @arr2;
	return grep( $arr1{$_}, @arr2 );
}

# for two arrays A and B, returns the array elements that only exist in array A and not in array B
sub _array_minus {
	my ($a, $b) = @_;
	my @a = @{$a};
	my @b = @{$b};
	my %a = map{$_=>1} @a;
	my %b = map{$_=>1} @b;
	return grep ( ! defined $b{$_}, @a );
}

# writes a tree description that spans provided taxa/clade to file
sub write_tree {
}

# Changes NCBI taxon identifiers in a given tree to taxon names
sub remap_to_name {
        my ($self, $tree) = @_;
        $tree->visit(sub{
                my $n = shift;
                if ( $n->is_terminal ) {
                        my $id = $n->get_name;
                        my $dbnode = $self->find_node($id);
                        $log->fatal("Could not find name for taxon id $id in database!") if not $dbnode;
                        my $name = $dbnode->taxon_name;
                        $name =~ s/_/\|/g;
                        $n->set_name( $name );
                }
                     });
        return $tree;
        
}

# Changes taxon names in a given tree to taxon identifiers
sub remap_to_ti { 
        my ($self, $tree) = @_;
        $tree->visit(sub{
                my $n = shift;
                if ( $n->is_terminal ) {
                        my $name = $n->get_name;
                        $name =~ s/_/ /g;                        
                        $name =~ s/\|/_/g;
                        $name =~ s/^'(.*)'$/$1/g;
                        $name =~ s/^"(.*)"$/$1/g;                        
                        my @nodes = $self->search_node({taxon_name=>$name})->all;
						$log->warn("found more than one database entry for taxon $name, using first entry.") if scalar (@nodes) > 1;						
						die "could not find database entry for taxon name $name " if scalar (@nodes) == 0;						
                        $n->set_name( $nodes[0]->ti );
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
		$counter++ if /^\s*tree\s/i; #*
	}	
	
	my $babs = int( $burnin * $counter );
	
	# create temporary outfile name
	my ( $outfh, $outfile ) = tempfile();
	close $outfh;
	
	# execute command
	my $tmpl = '%s -burnin %i %s %s 2> /dev/null';
	my $command = sprintf $tmpl, $config->TREEANNOTATOR_BIN, $babs, $infile, $outfile;
        $log->debug("running command $command");
        system($command) and die "Error building consensus: $?";
        my $newicktree = $self->parse_newick_from_nexus( $outfile );
        return $newicktree;
}

sub parse_newick_from_nexus {
        my ($self, $nexusfile) = @_;
        my $newick = "";
        open my $fh, '<', $nexusfile or die $!;
        my @lines = <$fh>;
        close $fh;
        my @rev = reverse @lines;
        if (scalar (@rev) > 1) { 
                $newick = $rev[1];
        }
        else {
                $log->warn("nexus file has less than two lines")
        }
        $newick =~ s/^tree.+\[\&R\]\s+//g;
        
        # get ids for taxon labels from 'Translate' section in nexus file
        my @sub;
        foreach (@lines) {
                push(@sub, split) if (/Translate/ .. /\;/);
        }    
        shift @sub;
        pop @sub;
        my %id_map = @sub;
        
        # remove trailing commas from idenifiers, if any
        for ( values %id_map ) { s/,$//g };
        
        # remove comments [things between square brackets] but keep the 'posterior' from comment
        #  and set the posterior as node name
        my @comments = ( $newick =~ m/(\[.+?\])/g );
  		
  		# iterate through comments and parse out posterior value, if given
        foreach my $c (@comments){
       	my @matches = ( $c =~ m/posterior=([0-9.]+)/g);
  			$log->warn("Found more than one posterior in comment tag") if scalar(@matches) > 1;
  			my $posterior = $matches[0] || "";
  			my $quoted = quotemeta $c;
  			
  			# substitute in newick tree string
  			$newick =~ s/$quoted/$posterior/g;
        }
		
        # substitute nexus taxon ids with real taxon labels
        # note that there is possible trouble if node names for posteriors (e.g. 1) 
        #   overlap with nexus identifier. However, BEAST seems to write all posteriors
        #   as proper decimals
		my $newicktree =  parse_tree(
                '-string'   => $newick,
                '-format' => 'newick',
                '-as_project' => 1,
        );

  	    $newicktree->visit( sub{
                my $n = shift;
                $n->set_name( $id_map{$n->get_name} ) if exists $id_map{$n->get_name};  
        });
        
        return $newicktree;
        
}

sub write_newick_tree {
        my ($self, $tree) = @_;
        my $str = "";
        $tree->visit_depth_first(
                '-pre_daughter'   => sub { $str .= '('             },     
                '-post_daughter'  => sub { $str .= ')'             },     
                '-in'             => sub { 
                	my $n = shift; 
                	$str .= $n->get_name;  
                	$str .= ":"  . $n->get_branch_length if $n->get_branch_length;},
                '-pre_sister'     => sub { $str.= ','             },     
            );
        $str .= ';';
        return $str;

}

# grafts a clade tree into a backbone tree, returns backbone
sub graft_tree {
	my ( $self, $backbone, $clade ) = @_;
    my $num_terminals = scalar @{ $backbone->get_terminals };
 
    # sometimes single quotes are added in the beast output when dealing with special characters
    #   removing them to match names in the backbone. Just to be sure, remove quotes also from backbone
    $clade->visit(sub{
                my $node = shift;
                my $name = $node->get_name;
                $name =~ s/\'//g;
                $node->set_name($name);
                      });
    $backbone->visit(sub{
                my $node = shift;
                my $name = $node->get_name;
                $name =~ s/\'//g;
                $node->set_name($name);
                         });
        
    my @names = map { $_->get_name } @{ $clade->get_terminals };
	$log->debug("Clade terminals : @names");
	
	# retrieve the tips from the clade tree that also occur in the backbone, i.e. the 
	# exemplars, and locate their MRCA on the backbone
	my @exemplars;
	for my $name ( @names ) {
		if ( my $e = $backbone->get_by_name($name) ) {
			$log->info("found exemplar $name in backbone: ".$e->get_name);
                        push @exemplars, $e;
		}
		else {
			$log->debug("$name is not in the backbone tree");
		}
	}
	$log->info("found ".scalar(@exemplars)." exemplars in backbone");
	my @copy = @exemplars; # XXX ???	
	my $bmrca = $backbone->get_mrca(\@copy);
	my $nodes_to_root = $bmrca->calc_nodes_to_root;
	$log->info("backbone MRCA distance from root : " . $nodes_to_root);
	if ( $bmrca->is_root ){
		$log->fatal("Something goes wrong here: MRCA of exemplar species " . join(',', map{$_->id} @exemplars) . " in backbone is the backbone root!");
	}
    # find the exemplars in the clade tree and find their MRCA
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
                $log->debug("adjusting branch length for ".$node->get_internal_name." to ".$bl * ( $bmrca_depth / $cmrca_depth ));
                      });
	$log->debug("re-scaled clade by $bmrca_depth / $cmrca_depth");
	
	# calculate the depth of the root in the clade, adjust backbone depth accordingly
	my $crd  = $clade->get_root->calc_max_path_to_tips;
	my $diff = $bmrca_depth - $crd;
        my $branch_length = $bmrca->get_branch_length || 0;
        $bmrca->set_branch_length( $branch_length + $diff );
		$log->debug("adjusted backbone MRCA depth by $bmrca_depth - $crd");
        # now graft!
        $log->debug("Backbone tree before grafting: \n " . $backbone->to_newick);
        $bmrca->clear();
        $clade->visit(
		sub {
			my $node = shift;
			my $name = $node->get_internal_name;
			if ( my $p = $node->get_parent ) {
				if ( $p->is_root ) {
					$log->info("grafted node with id $name onto backbone MRCA");
					$node->set_parent($bmrca);								
				}
			}
			if ( $node->is_root ) {
				$log->info("replacing root with id $name with backbone MRCA");
				##$node->set_parent($bmrca);								

			}
			else {
				$backbone->insert($node);
			}
		}
            );
       $log->debug("Backbone tree after grafting: \n " . $backbone->to_newick);
            
    my $num_terminals_after = scalar @{ $backbone->get_terminals };
	$log->info("Grafting changed number of taxa from " . $num_terminals . " to " . $num_terminals_after);
	return $backbone;
}

sub make_phylip_binary {
	my ( $self, $phylip, $binfilename, $parser, $work_dir) = @_;
	$log->info("going to make binary representation of $phylip => $binfilename");	
        $log->info("using parser $parser");	
	my ( $phylipvolume, $phylipdirectories, $phylipbase ) = File::Spec->splitpath( $phylip );
        my $curdir = getcwd;
	chdir $work_dir;       
        
        # check if filename exists in working directory, if not take the full path
        my $phylipfile = -e $phylipbase ? $phylipbase : $phylip;
        
        my @command = ( $parser, 
		'-m' => 'DNA', 
		'-s' => $phylipfile, 
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
		if ( /^\s?(\S+)/ ) {                        
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

sub extract_clades {
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
                                
                                my ($genus) = map { $_->{'genus'} } grep { $_->{'species'} eq $id || 
                                										   $_->{'subspecies'} eq $id || 
                                										   $_->{'varietas'} eq $id || 
                                										   $_->{'forma'} eq $id} @records;
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
                my @s = map { $mt->get_species_and_lower_for_taxon( 'genus' => $_, @records ) } @g;
                push @set, \@s if scalar(@s) > 2;
        }
        
        # one special case can be that all genera are monotypic. 
        # in this case, all leaves are returned in one set to
        # be inferred together.
        if ( !%monophyletic and !%paraphyletic ) {
        	$log->info("all genera in backbone are monotypic, all species will form a single clade ");
        	my @s = map { $mt->get_species_and_lower_for_taxon( 'genus' => $_, @records ) } keys(%all);
        	push @set, \@s;
        } 
        
        return @set;        
}


1;
