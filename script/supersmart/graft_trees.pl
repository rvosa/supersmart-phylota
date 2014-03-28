#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::PhyLoTA::Service::TreeService;


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
        my @cladetrees;
        
        # parse backbone tree
        my $backbone_tree = parse_tree(
                '-file'   => $backbone,
                '-format' => 'newick',
                '-as_project' => 1,
            );

        if ( $cladetree and -e $cladetree ) {
                push @cladetrees, $cladetree;
        } 
        # if no clade treefile given, iterate over work dir and look for directories named cladeXXX        
        else {
                opendir my $dh, $workdir or die $!;
                while( my $entry = readdir $dh ) {                        
                        if ( $entry =~ /clade\d+/ && -d "${workdir}/${entry}" ) {                                
                                # this should be a nexus file			
                                my $stem = "${workdir}/${entry}/${entry}";
                                my $file = "${stem}.nex";
                                if ( -e $file ) {
                                        push @cladetrees, $file;
                                }                                
                        }                        
                }       
        }
        # graft all clade trees onto backbone tree
        my $grafted = $backbone_tree;
        for my $treefile ( @cladetrees ) {
                my $consensus = $ts->consense_trees( '-infile' => $treefile );
                print "Grafted : ".ref($grafted)." Consensus : ".$consensus."\n";
                $grafted = $ts->graft_tree( $grafted, $consensus );                
        }
        
        # save final tree in newick format
        my $grafted_file = $outfile || "grafted.newick";
        if ( $workdir ) {
                $grafted_file = "${workdir}/${grafted_file}";
        }
        open my $outfh, '>', $grafted_file or die $!;
        print $outfh $grafted->to_newick;
        close $outfh;
}
else {
	$logger->fatal("need -workdir or -cladetree and -backbone arguments");
	exit 1;
}


