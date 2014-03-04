package Bio::Phylo::PhyLoTA::Service::TreeService;
use strict;
use warnings;
use File::Temp 'tempfile';
use Bio::Phylo::PhyLoTA::Config;

my $config = Bio::Phylo::PhyLoTA::Config->new;

# writes a tree description that spans provided taxa/clade to file
sub write_tree {

}

# it appears (obviously) that midpoint rooting is going to fail in some
# cases. How to automate this?
sub reroot_tree {
	
}

# -infile, -format, -outfile, -burnin
sub consense_trees {
	my ( $self, %args ) = @_;

}

# splices a clade tree into a backbone tree
sub graft_trees {

}

# builds a tree using the configured tree inference method
sub infer_backbone_tree {

}

sub infer_clade_tree {

}

1;
