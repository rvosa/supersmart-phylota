# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::Domain::BigTree;
use strict;
use warnings;
use Cwd;
use File::Copy;
use File::Basename;
use Bio::Phylo::PhyLoTA::Config;
use Bio::Phylo::Util::Logger;
use Bio::Phylo::Util::Exceptions 'throw';

my $conf = Bio::Phylo::PhyLoTA::Config->new;
my $log  = Bio::Phylo::Util::Logger->new;
my $cwd  = getcwd();

=head1 NAME

Bio::Phylo::PhyLoTA::Domain::BigTree - Big Tree

=head1 DESCRIPTION

Represents a Very Large Phylogeny. The idea, tentatively, is that there might be 
situations where commonly-used memory-based tree classes (e.g. 
L<Bio::Phylo::Forest::Tree> or L<Bio::Tree::TreeI>) are not scalable enough for the 
pipeline's needs, so this package would provide some kind of file-based way of interacting 
with big trees. In fact, this idea hasn't really come to fruition yet.

=head1 METHODS

=over

=item new

The constructor takes no arguments.

=cut

sub new {
    my $class = shift;
    my $self = bless {}, $class;
    return $self;
}

=item build_tree

Infers a tree, returns the location of the produced tree file. XXX: this functionality
should be moved to a subclass of L<Bio::Phylo::PhyLoTA::Service>.

=cut

sub build_tree {
    my $self = shift;
    my %args = @_;
    
    # these are required arguments
    my $seqfile  = $args{'-seq_file'}  or throw 'BadArgs' => "need -seq_file argument";
    my $treefile = $args{'-tree_file'} or throw 'BadArgs' => "need -tree_file argument";
    my $seqfile_basename = basename($seqfile);
    my $treefile_basename = basename($treefile);
    
    # these are normally derived from the config file
    my $work_dir   = $args{'-work_dir'}   || $conf->WORK_DIR;
    my $examl_bin  = $args{'-examl_bin'}  || $conf->EXAML_BIN;
    my $parser_bin = $args{'-parser_bin'} || $conf->PARSER_BIN;
    my $examl_args = $args{'-examl_args'} || $conf->EXAML_ARGS;
    my $mpirun_bin = $args{'-mpirun_bin'} || $conf->MPIRUN_BIN;
    my $nodes      = $args{'-nodes'}      || $conf->NODES;
    
    # create outfile from infile, if not provided
    my $outfile = $args{'-outfile'} || "${seqfile_basename}.dnd";
    
    # first create the binary representation of the seq file, if it doesn't exist yet
    if ( not -e "${work_dir}/${seqfile_basename}.binary" ) {
        $log->info("going to create binary compressed alignment file ${work_dir}/${seqfile_basename}.binary");
        chdir $work_dir;
        my $command = "$parser_bin -s $seqfile_basename -n $seqfile_basename -m DNA > /dev/null";
        $log->info("going to run '$command'");
        system( $command ) == 0 or throw 'BadArgs' => $?;
        chdir $cwd;
    }
    else {
        $log->info("binary compressed alignment file ${work_dir}/${seqfile_basename}.binary already existed");
    }
    
    # copy the starting tree to the working dir if it's not already there
    if ( not -e "${work_dir}/${treefile_basename}" ) {
        $log->info("copying $treefile to $work_dir");
        copy( $treefile, "${work_dir}/${treefile_basename}" );
    }
    
    # now run the tree search
    {
        chdir $work_dir;
        my $command = "$mpirun_bin -np $nodes $examl_bin $examl_args -s ${seqfile_basename}.binary -t $treefile_basename -n $outfile > /dev/null";
        $log->info("going to run '$command'");
        system( $command ) == 0 or throw 'BadArgs' => $?;
        chdir $cwd;
    }
    return "$work_dir/ExaML_result.$outfile";
}

=back

=cut


1;