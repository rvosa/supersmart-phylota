#!/bin/bash

# this shell script documents the steps of the pipeline in order. the idea
# is that there will be a smarter way of managing the workflow, e.g. using
# galaxy or taverna, so this is not the definitive pipeline. it's meant to
# show how things are supposed to fit together.

# project-specific variables
WORKDIR=examples/primates          # working directory for intermediate files
NAMELIST=$WORKDIR/names.txt        # input list of taxon names
FOSSILTABLE=$WORKDIR/fossils.tsv

# template for looking up variables in the config file
PRINTVAL="perl -MBio::Phylo::PhyLoTA::Config=printval -e printval"

# config vars, optionally change by setting env var with SUPERSMART prefix, e.g.
# $ SUPERSMART_WORK_DIR=/tmp
VERBOSE=`$PRINTVAL VERBOSITY`       # global verbosity level
PARSER=`$PRINTVAL PARSER_BIN`       # converts phylip to examl input files
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
MEGATREE=megatree.dnd
CHRONOGRAM=$WORKDIR/chronogram.dnd
TREEPLCONF=$WORKDIR/treePL.conf
FINALTREE=$WORKDIR/final.dnd

# creates a table where the first column has the input species
# names and subsequent columns have the species ID and higher taxon IDs
if [ ! -e $SPECIESTABLE ]; then
	$MPIRUN $PERLSCRIPT/mpi_write_taxa_table.pl -i $NAMELIST \
	$VERBOSE > $SPECIESTABLE
fi

# creates the NCBI common tree from an input species table
if [ ! -e $COMMONTREE ]; then
	$PERLSCRIPT/write_common_tree.pl --nodelabels -i $SPECIESTABLE \
	$VERBOSE > $COMMONTREE
fi

# create alignments from an input species table
if [ ! -e $ALIGNMENTLIST ]; then
	$MPIRUN $PERLSCRIPT/mpi_write_alignments.pl -i $SPECIESTABLE \
	$VERBOSE -w $WORKDIR > $ALIGNMENTLIST
fi

# merge alignments by orthology 
if [ ! -e $MERGEDLIST ]; then
	$PERLSCRIPT/merge_alignments.pl -l $ALIGNMENTLIST -w $WORKDIR \
	$VERBOSE > $MERGEDLIST
fi

# create supermatrix for backbone taxa
if [ ! -e $SUPERMATRIX ]; then
	$PERLSCRIPT/pick_exemplars.pl -l $MERGEDLIST -t $SPECIESTABLE $VERBOSE > $SUPERMATRIX
fi

# infer backbone tree
if [ ! -e "$WORKDIR/$MEGATREE" ]; then
	$PERLSCRIPT/infer_backbone.pl -o "$WORKDIR/$MEGATREE" -w $WORKDIR -c $COMMONTREE \
	-s $SUPERMATRIX $VERBOSE
fi

# create treePL config file
if [ ! -e $TREEPLCONF ]; then
	NUMSITES=`head -1 $SUPERMATRIX | cut -f 2 -d ' '`
	$PERLSCRIPT/write_treepl_config.pl -f $FOSSILTABLE -r "$WORKDIR/$MEGATREE" \
	-s $TREEPLSMOOTH -w $CHRONOGRAM $VERBOSE -n $NUMSITES > $TREEPLCONF
fi

# run treePL
if [ ! -e $CHRONOGRAM ]; then
    $TREEPLBIN $TREEPLCONF
fi

# decompose backbone tree
if [ ! -e $WORKDIR/clade* ]; then
    $PERLSCRIPT/decompose_backbone.pl -t $SPECIESTABLE -l $MERGEDLIST -b "$WORKDIR/$MEGATREE" -w $WORKDIR $VERBOSE
fi

# merge clade alignments
clade_count=`ls -d $WORKDIR/clade*/ 2>/dev/null | wc -w`
file_count=`ls -d $WORKDIR/clade*/*.xml 2>/dev/null | wc -w`
if [ $file_count -le $clade_count ]; then    
    echo "HIER"
    $PERLSCRIPT/merge_clade_alignments.pl -w $WORKDIR $VERBOSE
fi

# infer tree for each clade
file_count=`ls -f $WORKDIR/clade*/*.nex 2>/dev/null | wc -w`
if [ $file_count -le $clade_count ]; then
    echo "HIER2"
    $PERLSCRIPT/infer_clade.pl -w $WORKDIR -ngens 30000000 -sfreq 300000 -lfreq 300000 $VERBOSE
fi

# graft clade trees onto backbone
if [ ! -e $FINALTREE ]; then
    $PERLSCRIPT/graft_trees.pl -workdir $WORKDIR -backbone $NAMED_CHRONOGRAM -outfile $FINALTREE $VERBOSE
fi
