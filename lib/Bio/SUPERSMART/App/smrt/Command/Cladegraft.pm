package Bio::SUPERSMART::App::smrt::Command::Cladegraft;

use strict;
use warnings;

use Bio::Phylo::PhyLoTA::Service::TreeService;

use Bio::Phylo::IO 'parse_tree';

use base 'Bio::SUPERSMART::App::smrt::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: grafts the inferred clade trees on the backbone chronogram

=head1 NAME

Cladegraft.pm - Grafts the inferred clade trees on the backbone chronogram.

=head1 SYNOPSYS

smrt cladegraft [-h ] [-v ] [-w <dir>] -b <file> [-o <file>] [-c <file>]

=head1 DESCRIPTION

Combines a backbone tree of representative genera with one or more clade trees which have been
infered independently. Given a directory as argument 'workdir', traverses it, looks for subdirectories
and files that match the pattern clade*.nex. These must be NEXUS files. 
Given a single NEXUS file file as 'cladetree' argument, grafts this tree onto the backbone.
The resulting tree is exported in the NEWICK format.

=cut

sub options {
	my ($self, $opt, $args) = @_;
	my $outfile_default = "final.dnd";
	return (                               	
  		["backbone|b=s", "backbone tree as produced by 'smrt bbinfer'", { arg => "file", mandatory => 1}],
		["outfile|o=s", "name of the output tree file (newick format) defaults to $outfile_default", { default=> $outfile_default, arg => "file"}],    	    
		["cladetree|c=s", "name of tree file for single clade (newick format) if only a single cladetree should be grafted", { arg => "file"}],    	    		
    );
}

sub validate {
	my ($self, $opt, $args) = @_;		

	if ( ! $opt->backbone) {
		$self->usage_error("need backbone argument!");
	}

	if ( my  $file = $opt->backbone ){
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
	if ( my  $file = $opt->cladetree ) {
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);							
	}

	# check if there exist 'clade folders' in working directory
	if (! $opt->cladetree){
		(my $workdir = $opt->workdir) =~ s/\/$//g;
		if ( ! -d glob "$workdir/clade*" ){			
        	$self->usage_error("no cladetree argument given and no clade folders found in working directory $workdir");
		}
	}
}

sub run{
	my ($self, $opt, $args) = @_;		
	
	# collect command-line arguments
	my $backbone = $opt->backbone;
	my $cladetree = $opt->cladetree;
	my $outfile = $self->outfile;
	(my $workdir = $opt->workdir) =~ s/\/$//g;
 
	# instantiate helper objects
	my $ts = Bio::Phylo::PhyLoTA::Service::TreeService->new;
	my $logger = $self->logger;
	
	# parse backbone tree
	my $backbone_tree = parse_tree(
                '-file'   => $backbone,
                '-format' => 'newick',
                '-as_project' => 1,
            );
		
	# filenames for nexus trees written by BEAST
    
    
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
    open my $outfh, '>', $outfile or die $!;
    print $outfh $ts->write_newick_tree($grafted);
    close $outfh;	

	$logger->info("DONE, results written to $outfile");
}

1;