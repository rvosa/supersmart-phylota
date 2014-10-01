#!/bin/bash

# this shell script documents the steps of the pipeline in order. the idea
# is that there will be a smarter way of managing the workflow, e.g. using
# galaxy or taverna, so this is not the definitive pipeline. it's meant to
# show how things are supposed to fit together.

# project-specific variables
WORKDIR=examples/fishes         # working directory for intermediate files
NAMELIST=$WORKDIR/names.txt 	   # input list of taxon names
FOSSILTABLE=$WORKDIR/fossils.tsv

# template for looking up variables in the config file
PRINTVAL="perl -MBio::Phylo::PhyLoTA::Config=printval -e printval"

# config vars, optionally change by setting env var with SUPERSMART prefix, e.g.
# $ SUPERSMART_WORK_DIR=/tmp
VERBOSE=`$PRINTVAL VERBOSITY`       # global verbosity level
PERL=`$PRINTVAL PERL_BIN`           # the perl interpreter
MPIBIN=`$PRINTVAL MPIRUN_BIN`       # MPI job dispatcher, 'mpirun' or 'mpiexec'
NODES=`$PRINTVAL NODES`             # number of nodes and/or threads
MPIRUN="$MPIBIN -np $NODES"         # invocation for parallel jobs
EXAMLBIN=`$PRINTVAL EXAML_BIN`      # location of examl
EXAMLARGS=`$PRINTVAL EXAML_ARGS`    # command line arguments for examl
EXAML="$EXAMLBIN $EXAMLARGS"        # invocation of examl
TREEPLBIN=`$PRINTVAL TREEPL_BIN`    # location of treePL
TREEPLSMOOTH=`$PRINTVAL TREEPL_SMOOTH`

# shorthand for running the perl scripts without having to specify the added
# search path (the lib folder) and the location of the scripts
PERLSCRIPT="$PERL $SUPERSMART_HOME/script/supersmart"

# names for intermediate files within the working directory
SPECIESTABLE=$WORKDIR/species.tsv
COMMONTREE=$WORKDIR/common.dnd
ALIGNMENTLIST=$WORKDIR/aligned.txt
MERGEDLIST=$WORKDIR/merged.txt
SUPERMATRIX=$WORKDIR/supermatrix.phy
MEGATREE=$WORKDIR/megatree.dnd
REROOTED_MEGATREE=$WORKDIR/rerooted_megatree.dnd
BACKBONE=$WORKDIR/backbone.dnd
TREEPLCONF=$WORKDIR/treePL.conf
FINALTREE=$WORKDIR/final.dnd

# creates a table where the first column has the input species
# names and subsequent columns have the species ID and higher taxon IDs
if [ ! -e $SPECIESTABLE ]; then
	$MPIRUN $PERLSCRIPT/parallel_write_taxa_table.pl -i $NAMELIST \
	$VERBOSE > $SPECIESTABLE
fi

# creates the NCBI common tree from an input species table
if [ ! -e $COMMONTREE ]; then
	$PERLSCRIPT/write_common_tree.pl --nodelabels -i $SPECIESTABLE \
	$VERBOSE > $COMMONTREE
fi

# create alignments from an input species table
if [ ! -e $ALIGNMENTLIST ]; then
	$MPIRUN $PERLSCRIPT/parallel_write_alignments.pl -i $SPECIESTABLE \
	$VERBOSE -w $WORKDIR > $ALIGNMENTLIST
fi

# merge alignments by orthology 
if [ ! -e $MERGEDLIST ]; then
	$MPIRUN $PERLSCRIPT/parallel_merge_alignments.pl -l $ALIGNMENTLIST -w $WORKDIR \
	$VERBOSE > $MERGEDLIST
fi

# create supermatrix for backbone taxa
if [ ! -e $SUPERMATRIX ]; then
	$PERLSCRIPT/pick_exemplars.pl -l $MERGEDLIST -t $SPECIESTABLE $VERBOSE > $SUPERMATRIX
fi


# infer backbone tree
if [ ! -e $MEGATREE ]; then
	$PERLSCRIPT/infer_backbone.pl -o $MEGATREE -w $WORKDIR -c $COMMONTREE -s $SUPERMATRIX $VERBOSE
fi


# reroot backbone tree
if [ ! -e $REROOTED_MEGATREE ]; then
	$PERLSCRIPT/reroot_backbone.pl -w $WORKDIR -b $MEGATREE -t $SPECIESTABLE $VERBOSE > $REROOTED_MEGATREE
fi

# remap to taxon names
if [ ! -e $BACKBONE ]; then
	$PERLSCRIPT/remap.pl $REROOTED_MEGATREE > $BACKBONE
fi

# decompose backbone tree
count=`ls -d $WORKDIR/clade* 2>/dev/null | wc -w`
if [ $count -eq 0 ]; then
    $MPIRUN $PERLSCRIPT/parallel_decompose_backbone.pl -t $SPECIESTABLE -l $MERGEDLIST -b $REROOTED_MEGATREE -w $WORKDIR $VERBOSE
fi

# merge clade alignments
clade_count=`ls -d $WORKDIR/clade*/ 2>/dev/null | wc -w`
file_count=`ls -d $WORKDIR/clade*/*.xml 2>/dev/null | wc -w`

if [ $file_count -lt $clade_count ]; then    
        $PERLSCRIPT/merge_clade_alignments.pl -w $WORKDIR $VERBOSE
fi

# infer tree for each clade
file_count=`ls -f $WORKDIR/clade*/*.nex 2>/dev/null | wc -w`
if [ $file_count -lt $clade_count ]; then
  	 $MPIRUN $PERLSCRIPT/parallel_infer_clades.pl -w $WORKDIR -ngens 3000000 -sfreq 30000 -lfreq 30000 $VERBOSE
fi

# graft clade trees onto backbone
if [ ! -e $FINALTREE ]; then
    $PERLSCRIPT/graft_trees.pl -workdir $WORKDIR -backbone $BACKBONE -outfile $FINALTREE $VERBOSE
fi
