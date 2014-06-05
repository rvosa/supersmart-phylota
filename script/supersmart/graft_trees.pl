#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Service::TreeService;

=head1 NAME

graft_trees.pl - fuses backbone tree with one or more clade trees

=head1 SYNOPSYS
 
 $ perl graft_trees.pl -workdir <directory> -backbone <file> -cladetree <file> \
    -outfile <file> [--verbose]
    
=head1 DESCRIPTION

Combines a backbone tree of representative genera with one or more clade trees which have been
infered independently. Given a C<directory>  as argument C<-workdir>, traverses it, looks for subdirectories
and files that match the pattern C<clade\d+/clade\d+\.nex>. These must be NEXUS files. 
Given a single NEXUS C<file> file as C<-cladetree> argument, grafts this tree onto the backbone.
The resulting tree is exported in the NEWICK format.

=over

=item outfile

Optional file name to which the final tree is written. Defaults to "grafted.dnd".

=item backbone

Name of the file containing the backbone tree in NEWICK format.

=cut 

# process command line objects
my $verbosity = WARN;
my $workdir;
my $backbone;
my $cladetree;
my $outfile;

GetOptions(
	'workdir=s'  => \$workdir,
	'backbone=s' => \$backbone,
        'cladetree=s'=> \$cladetree,
        'outfile=s'  => \$outfile,
        'verbose+'   => \$verbosity,        
        
);

# instantiate helper objects
my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
my $logger = Bio::Phylo::Util::Logger->new( 
	'-level' => $verbosity, 
	'-class' => [ 
		'main',
		'Bio::Phylo::Phylota::Service::Treeservice',
	],
);

if ( (($workdir and -d $workdir) or ($cladetree and -e $cladetree)) and $backbone and -e $backbone ) {
        # filenames for clade trees
        my %cladetrees;
        
        # parse backbone tree
        my $backbone_tree = parse_tree(
                '-file'   => $backbone,
                '-format' => 'newick',
                '-as_project' => 1,
                #'-ignore_comments' =>0,
            );

        #if ( $cladetree and -e $cladetree ) {
        #        push @cladetrees, $cladetree;
        #} 
        # if no clade treefile given, iterate over work dir and look for directories named cladeXXX        
        my $grafted = $backbone_tree;
        
        opendir my $dh, $workdir or die $!;
        while( my $entry = readdir $dh ) {                        
                if ( $entry =~ /clade\d+/ && -d "${workdir}/${entry}" ) {                                
                        # this should be a nexus file			
                        my $stem = "${workdir}/${entry}/${entry}";
                        my $file = "${stem}.nex";
                        if ( -e $file ) {
                                $logger->info( "Processing clade $entry" );

                                my $consensus = $ts->consense_trees( '-infile' => $file ); 

                                # save consensus tree
                                open my $fh, '>', $stem.".dnd" or die $!;
                                print $fh $ts->write_newick_tree($consensus);
                                close $fh;
                                $grafted = $ts->graft_tree( $grafted, $consensus );                
                        }                               
                }                        
        }       

        # save final tree in newick format
        my $grafted_file = $outfile || "grafted.dnd";
        open my $outfh, '>', $grafted_file or die $!;
        print $outfh $ts->write_newick_tree($grafted);
        close $outfh;
}
else {
	$logger->fatal("need -workdir or -cladetree and -backbone arguments");
	exit 1;
}


