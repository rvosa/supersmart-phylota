package Bio::Phylo::PhyLoTA::Service::TreeService;
use strict;
use warnings;
use Cwd;
use File::Temp 'tempfile';
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Factory;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::PhyLoTA::Service;
use base 'Bio::Phylo::PhyLoTA::Service';
use Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa;
use List::MoreUtils 'uniq';
use List::Util 'min';

my $config = Bio::Phylo::PhyLoTA::Config->new;
my $log = Bio::Phylo::PhyLoTA::Service->logger->new;	
my $fac = Bio::Phylo::Factory->new;

=head1 NAME

Bio::Phylo::PhyLoTA::Service::TreeService - Operations on phylogenetic Trees

=head1 DESCRIPTION

Provides functionality for handling trees, e.g. decomposing, rerooting and grafting
clade trees onto a backbone tree.


=over


=item reroot_tree

Reroots a backbone tree. All possible rerootings are evaluated an the tree is returned that
minimizes the paraphyly with respect to the exemplar species of all genera.

=cut

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
	my $newi          = $tree->to_newick(nodelabels=>1);

	# iterate over all possible rerootings
	for my $i (0 .. $num_internals - 1 ) {

		# get fresh tree, so the node order etc. won't be messed up
		my $current_tree = parse(
			'-string' => $newi,
			'-format' => 'newick',
		)->first;

		my @internals = @{ $current_tree->get_internals };

		# reroot the tree
		if ( my $node = $internals[$i] ) {
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
		$log->warn ("Found more than one optimal rerooted tree. Returning the first one");
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


=item remove_internal_names 

Removes those labels of internal and nodes in a tree, that are not numeric, e.g. "root" of
"e1" which may have been introduced by rerooting or resolving the tree

=cut

sub remove_internal_names {
	my ($self, $tree) = @_;	
    
    for my $n ( @{$tree->get_internals} ) {
		if ($n->get_name =~ m/[a-zA-Z]/){
			$n->set_name('');
		}
    }
    my $n = $tree->get_root;
    if ($n->get_name =~ m/[a-zA-Z]/){
		$n->set_name('');
	}
	return $tree;
}

=item remap_to_name

Given an object of class L<Bio::Phylo::Forest::Tree>, 
changes the names of the terminal nodes from NCBI taxonomy identifiers to their respective taxon names. 

=cut


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

=item remap_to_ti

Given an object of class L<Bio::Phylo::Forest::Tree>, 
changes the names of all terminal nodes from taxon names to their respective identifiers 
as given in the NCBI taxonomy database. 

=cut

sub remap_to_ti { 
        my ($self, $tree) = @_;

        $tree->visit(sub{
                my $n = shift;
                if ( $n->is_terminal and my $name = $n->get_name ) {      
                        $name =~ s/_/ /g;                        
                        $name =~ s/\|/_/g;
                        $name =~ s/^'(.*)'$/$1/g;
                        $name =~ s/^"(.*)"$/$1/g;                        
						my $dbnode = $self->find_node({taxon_name=>$name});                     
						die "could not find database entry for taxon name $name " if not $dbnode;						
         				my $ti = $dbnode->ti;
                        $n->set_name( $ti );						
                        $log->debug("Remapped name $name to $ti ");
                }
                     });
        return $tree;       
}

=item consense_trees

Given a file with multiple trees (in newick format), 
creates a consensus tree using the tool specified by 
"TREEANOTATOR_BIN" in the configuration file.
Returns a single newick tree as a string.

=cut

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

=item parse_newick_from_nexus

Given the file name of a phylogenetic tree in Nexus format, creates a tree in newick representation which is
returned as a string. If information about the posterior value of a node is given (as present in the Nexus output of *BEAST),
the posterior value is assigned to the respective inner nodes in the newick tree.



=cut 

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
		my $newicktree =  parse(
                '-string'   => $newick,
                '-format' => 'newick',
        )->first;

  	    $newicktree->visit( sub{
                my $n = shift;
                $n->set_name( $id_map{$n->get_name} ) if exists $id_map{$n->get_name};  
        });
        
        return $newicktree;
        
}

=item graft_tree

Grafts a clade tree into a backbone tree, returns the altered backbone.

=cut

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
    $log->debug("Backbone tree before grafting: \n " . $backbone->to_newick);
    $log->debug("Depth of exemplar mrca in backbone : $bmrca_depth");
    $log->debug("Depth of exemplar mrca in clade : $cmrca_depth");
        
    $clade->visit(sub{
		my $node = shift;
		my $bl = $node->get_branch_length || 0; # zero for root
        if ( $bmrca_depth > 0 and $cmrca_depth > 0 ){  		
       		$node->set_branch_length( $bl * ( $bmrca_depth / $cmrca_depth ) );
        	$log->debug("adjusting branch length for ".$node->get_internal_name." to ".$bl * ( $bmrca_depth / $cmrca_depth ));
        } 
        else {
        	$log->warn("mrca of clade or backbone has depth zero!");
        	}
    });
		
	#$log->debug("re-scaled clade by $bmrca_depth / $cmrca_depth");
	
	# calculate the depth of the root in the clade, adjust backbone depth accordingly
	my $crd  = $clade->get_root->calc_max_path_to_tips;
	my $diff = $bmrca_depth - $crd;
        my $branch_length = $bmrca->get_branch_length || 0;
        $bmrca->set_branch_length( $branch_length + $diff );
		$log->debug("adjusted backbone MRCA depth by $bmrca_depth - $crd");
        # now graft!
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

=item make_phylip_binary

Given the location of a file in phylip format, creates a binary representation of
the phylip file and writes it to the specified file. 

=cut

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

=item make_phylip_from_matrix

Given an object of class L<Bio::Phylo::Matrices::Datum>, writes the data into a phylip
file to the specified file name.

=cut

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


=item read_tipnames

Reads the supermatrix (phylip format) and returns the tip names from it

=cut

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

=item extract_clades

Given a backbone tree and a list with the taxa classifications, 
decomposes the backbone tree into single clades, represented by lists of 
taxon ids. The decomposition identifies monophyletic genera which will form individual clades.
If genera are paraphyletic in the tree (i.e. the the mrca of two exemplar species in one genus
is not the direct parent of the two exemplars), we traverse up the tree and include all species 
from sister clades until the paraphyly is resolved.

=cut

sub extract_clades {
	my ($self, $tree, @records) = @_;
	my %genera;
    my $mt = 'Bio::Phylo::PhyLoTA::Domain::MarkersAndTaxa';
    
    # first get all taxa that will be included in the clades
    my @valid_ranks = ("species", "subspecies", "varietas", "forma");
    my $level = $mt->get_highest_informative_level(@records);
    my @highest_taxa = $mt->get_distinct_taxa($level, @records);
    my %all_taxa = map{$_=>1 } $mt->query_taxa_table(\@highest_taxa, \@valid_ranks, @records);
    	    
	# get exemplars for genera
	my %terminal_ids = map {$_=>1} map{$_->get_name} @{ $tree->get_terminals };
    
	foreach my $id ( keys(%terminal_ids) ) {
		my ($genus) = $mt->query_taxa_table($id, "genus", @records);
		push @{$genera{$genus}->{'exemplars'}}, $id;			
	}
	    
	# get mrca and distance from root for all exemplars in each genus
	foreach my $genus ( keys %genera ) {	
		my %exemplar_ids = map {$_=>1} @{$genera{$genus}->{'exemplars'}};		
		my @nodes = grep { exists( $exemplar_ids{$_->get_name}) }  @{$tree->get_terminals};
		my $mrca = $tree->get_mrca( \@nodes );
		
		$genera{$genus}->{'mrca'} = $mrca;
		my $dist = $mrca->calc_nodes_to_root;
		$genera{$genus}->{'dist_from_root'} = $dist;			
	}
	# sort genera by depth of exemplar mrca
	my @sorted_genera = sort { $genera{$a}->{'dist_from_root'} <=> $genera{$b}->{'dist_from_root'} } keys (%genera);
	
	# make sets of leaves that are in the subtree spanned from the respective mrca
	my @sets;
	for my $genus ( @sorted_genera ) {
		my $mrca = $genera{$genus}->{'mrca'};
		my $dist = $genera{$genus}->{'dist_from_root'};

		# replace mrca with parent node to include monotypic genera in sister clade
		#if ($mrca->is_terminal) {
		#	$mrca = @{$mrca->get_ancestors}[0];
		#}

		# get genera within the subtree:
		my @mrca_terminal_ids = map{$_->get_name} @{$mrca->get_terminals};
		my @mrca_genera =  $mt->query_taxa_table(\@mrca_terminal_ids, 'genus', @records);
		
				
		# get all species that are contained in the genera of the mrca 
		my @s = $mt->query_taxa_table(\@mrca_genera, \@valid_ranks, @records);
		# only include taxa in a set if they are not already part of a clade spanned higher in the hierarchy 
		my @remaining_terminals = keys(%all_taxa);

		my @clade = _intersect(\@remaining_terminals, \@s);
		delete @all_taxa{@clade}; 
		
		push @sets, \@clade;
	}
	
	# All sets with only one taxon represent monotypic genera. Since we have resolved 
	#  paraphyly at this point, all monotypic genera nested in other clades are accounted for;
	#  therefore it is safe to remove all non-nested monotypic genera.
	
	# There is one special case, however: All genera are monotypic, so all species will be
	#  put together into one clade!
	my @num_elems = uniq map{scalar(@$_)} @sets;
	if ( scalar(@num_elems)==1 and $num_elems[0] == 1 ) {
		@sets = [ map{$_} @sets ];
	}
	else {
		# remove monotypic genera from clades and sets that have <3 species and thus cannot be resolved
		@sets = grep {scalar(@$_)>2} @sets;
	}
	$log->info("Extracted " . scalar(@sets) . " clades");
	return @sets;	
}

=back 

=cut

1;
