#!/bin/bash

# template for looking up variables in the config file
PRINTVAL="perl -Ilib -MBio::Phylo::PhyLoTA::Config=printval -e printval"

# config vars, optionally change by setting env var with SUPERSMART prefix, e.g.
# $ export SUPERSMART_WORK_DIR=/tmp
WORKDIR=`$PRINTVAL WORK_DIR`
NAMELIST=`$PRINTVAL NAME_LIST_FILE`
VERBOSE=`$PRINTVAL VERBOSITY`
PARSER=`$PRINTVAL PARSER_BIN`
MYSQL=`$PRINTVAL MYSQL_BIN`
PERL=`$PRINTVAL PERL_BIN`
CURL=`$PRINTVAL CURL_BIN`
WGET=`$PRINTVAL WGET_BIN`
MPIBIN=`$PRINTVAL MPIRUN_BIN`
NODES=`$PRINTVAL NODES`
MPIRUN="$MPIBIN $NODES"
EXAMLBIN=`$PRINTVAL EXAML_BIN`
EXAMLARGS=`$PRINTVAL EXAML_ARGS`
EXAML="$EXAMLBIN $EXAMLARGS"
GUNZIP=`$PRINTVAL GUNZIP_BIN`

# shorthand for running the perl scripts without having to specify the added
# search path (the lib folder) and the location of the scripts
PERLSCRIPT="$PERL -Ilib script/supersmart"

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

# creates a table where the first column has the input species
# names and subsequent columns have the species ID and higher taxon IDs
if [ ! -e $SPECIESTABLE ]; then
	$MPIRUN $PERLSCRIPT/mpi_write_taxa_table.pl -i $NAMELIST $VERBOSE > $SPECIESTABLE
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
	$PERLSCRIPT/join_alignments.pl -l $MERGEDLIST $VERBOSE --nexus > $SUPERMATRIX
fi

# compresses the supermatrix into binary format
if [ ! -e "$WORKDIR/${BINARYMATRIX}.binary" ]; then
	cd $WORKDIR/ && $PARSER -s supermatrix.phy -n supermatrix-bin -m DNA && cd -
fi

# create input tree for examl
if [ ! -e $USERTREE ]; then
	$PERLSCRIPT/write_constraint_tree.pl -t $COMMONTREE -s $SUPERMATRIX \
	$VERBOSE > $USERTREE
fi

# run examl
if [ ! -e "$WORKDIR/$MEGATREE" ]; then
	cd $WORKDIR/ && $EXAML -s "${BINARYMATRIX}.binary" -t $USERTREE -n $MEGATREE && cd -
fi

