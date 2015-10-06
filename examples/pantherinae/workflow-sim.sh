#!/bin/bash

# This script uses simulations to replicate a dataset from a pervious
# SUPERSMART run. Replication includes making a synthetic tree 
# resembling the final tree from the analysis, an appropriate taxa
# file and simulated sequence alignment files.
# Subsequently, the supersmart pipeline is ran again 
# on the synthetic dataset. The tree infered from this tree inference
# can then directly be compared to the 'original' replicated tree. Thereby,
# the performance of the SUPERSMART method can be validated. 
 
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

# Replicate the dataset (final tree, taxa table and alignments)
smrt-utils replicate -t final.nex -f nexus -a aligned.txt -l replicate.log -v

# make plot with alignment summary
Rscript ../plot-aln-summary.R aligned-replicated.txt

# insert simulated sequences and possible artificial taxa into the database
# prefix for sequence accessions is the simulation directory
smrt-utils dbinsert -s aligned-replicated.txt -t taxa-replicated.tsv -p $SIMDIR -v

# Rerun the supersmart pipeline in the directory of the replicated dataset
smrt orthologize -i aligned-smrt-inserted.txt
smrt bbmerge -t taxa-replicated.tsv -a merged.txt
smrt bbinfer -i exabayes -s supermatrix.phy
smrt bbreroot -b backbone.dnd -t taxa-replicated.tsv
smrt bbcalibrate -t backbone-rerooted.dnd -f fossils.tsv
smrt consense -i chronogram.dnd
smrt bbdecompose -b consensus.nex -a aligned-smrt-inserted.txt -t taxa-replicated.tsv
smrt clademerge --enrich
smrt cladeinfer --ngens=15_000_000 --sfreq=1000 --lfreq=1000
smrt cladegraft

# clean database from artificial sequences
sqlite3 $SUPERSMART_HOME/data/phylota.sqlite 'delete from seqs where acc like "$SIMDIR%"'
sqlite3 $SUPERSMART_HOME/data/phylota.sqlite 'delete from nodes_194 where common_name like "$SIMDIR%"'
