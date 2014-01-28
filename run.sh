#!/bin/bash

# trhis shell script documents the steps of the pipeline in order. the idea
# is that there will be a smarter way of managing the workflow, e.g. using
# galaxy or taverna, so this is not the definitive pipeline. it's meant to
# show how things are supposed to fit together.

# template for looking up variables in the config file
PRINTVAL="perl -MBio::Phylo::PhyLoTA::Config=printval -e printval"

# config vars, optionally change by setting env var with SUPERSMART prefix, e.g.
# $ SUPERSMART_WORK_DIR=/tmp
WORKDIR=`$PRINTVAL WORK_DIR`        # working directory for intermediate files
NAMELIST=`$PRINTVAL NAME_LIST_FILE` # input list of taxon names
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
FOSSILTABLE=`$PRINTVAL FOSSIL_TABLE_FILE`
TREEPLSMOOTH=`$PRINTVAL TREEPL_SMOOTH`

# shorthand for running the perl scripts without having to specify the added
# search path (the lib folder) and the location of the scripts
PERLSCRIPT="$PERL $PHYLOTA_HOME/script/supersmart"

# names for intermediate files within the working directory
SPECIESTABLE=$WORKDIR/species.tsv
COMMONTREE=$WORKDIR/common.dnd
USERTREE=$WORKDIR/user.dnd
ALIGNMENTLIST=$WORKDIR/aligned.txt
ALIGNMENTSTEM=$WORKDIR/aligned.fa.
MERGEDLIST=$WORKDIR/merged.txt
MERGEDSTEM=$WORKDIR/merged.fa.
SUPERMATRIX=$WORKDIR/supermatrix.phy
CHUNKTABLE=$WORKDIR/chunks.tsv
CHUNKLIST=$WORKDIR/chunks.txt
BINARYMATRIX=supermatrix-bin
MEGATREE=megatree.dnd
CHRONOGRAM=$WORKDIR/chronogram.dnd
TREEPLCONF=$WORKDIR/treePL.conf

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
	$MPIRUN $PERLSCRIPT/mpi_merge_alignments.pl -l $ALIGNMENTLIST \
	$VERBOSE -s $MERGEDSTEM > $MERGEDLIST
fi

# join alignments by taxon
if [ ! -e $SUPERMATRIX ]; then
	$PERLSCRIPT/join_alignments.pl -l $MERGEDLIST $VERBOSE > $SUPERMATRIX
fi

# compresses the supermatrix into binary format
if [ ! -e "$WORKDIR/${BINARYMATRIX}.binary" ]; then
	cd $WORKDIR
	$PARSER -s supermatrix.phy -n $BINARYMATRIX -m DNA
	cd -
fi

# create input tree for examl
if [ ! -e $USERTREE ]; then
	$PERLSCRIPT/write_constraint_tree.pl -t $COMMONTREE -s $SUPERMATRIX \
	$VERBOSE > $USERTREE
fi

# run examl
if [ ! -e "$WORKDIR/$MEGATREE" ]; then
	cd $WORKDIR
	$EXAML -s "${BINARYMATRIX}.binary" -t $USERTREE -n $MEGATREE
	cd -
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

