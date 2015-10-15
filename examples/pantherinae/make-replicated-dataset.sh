#!/bin/bash

# This script uses simulations to replicate a dataset from a previous
# SUPERSMART run. Replication includes making a synthetic tree 
# resembling the final tree from the analysis, an appropriate taxa
# file and simulated sequence alignment files.
 
# Note that running this script will take considerable time and resources,
# so it is best to run it on a cluster computer.

# Note that this script must be called from the working directory of a 
# previous (completed!) SUPERSMART run.

# A directory name to store all simulated data and analysis results can be provided
# as an argument, if not given it defaults to 'simulations/' within the 
# current workind directory

# parse command-line arguments for simulation directory name
if [ -z "$1" ]
	then
	SIMDIR='simulations'
else
	SIMDIR=$1
fi

# make directory to store all files for simulation analysis
mkdir $SIMDIR

# copy files from previous SUPERSMART run necessary to replicate tree and data
cp final.nex $SIMDIR
cp aligned.txt $SIMDIR
cp fossils.tsv $SIMDIR
cd $SIMDIR

# parameters for SUPERSMART run

# Replicate the dataset (final tree, taxa table and alignments)
smrt-utils replicate -t final.nex -f nexus -a aligned.txt -l replicate.log -v

# insert simulated sequences and possible artificial taxa into the database
# prefix for sequence accessions is the simulation directory
smrt-utils dbinsert -s aligned-replicated.txt -t taxa-replicated.tsv -p $SIMDIR -v

# make file with alignment statistics
smrt-utils  alnstats -a aligned-smrt-inserted.txt

