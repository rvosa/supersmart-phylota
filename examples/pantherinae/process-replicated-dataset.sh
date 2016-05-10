#!/bin/bash

# This script runs a SUPERSMART analysis on a simulated dataset.
# The directory in which this dataset is located can be given as
# argument. Filenames for alignment lists etc are the same as
# produced by the script make-replicated-dataset.sh

# Directory can be given as argument, defaults to current directory
if [ -z "$1" ]
	then
	SIMDIR=$PWD
else
	SIMDIR=$1
fi

cd $SIMDIR

# Run the supersmart pipeline in the directory of the replicated dataset
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

# Make supermatrix with all taxa for comparison with other tools
# For this, we set the backbone parameters to the clade parameters to
# make a large dataset
export SUPERSMART_BACKBONE_MAX_DISTANCE=$SUPERSMART_CLADE_MAX_DISTANCE
export SUPERSMART_BACKBONE_MIN_COVERAGE=$SUPERSMART_CLADE_TAXON_MIN_MARKERS
export SUPERSMART_BACKBONE_MAX_COVERAGE=$SUPERSMART_CLADE_MAX_MARKERS;

# have to orthologize again because maximum distance in alignment changed
rm cluster*.fa
smrt orthologize -i aligned-smrt-inserted.txt

# make one matrix for Mr. Bayes and one in phylip format
smrt bbmerge -t taxa-replicated.tsv -a merged.txt -e -1 -o supermatrix-all.phy
smrt bbmerge -t taxa-replicated.tsv -a merged.txt -e -1 -o supermatrix-all.nex -f mrbayes

# Clean database from artificial sequences
sqlite3 $SUPERSMART_HOME/data/phylota.sqlite 'delete from seqs where acc like "$SIMDIR%"'
sqlite3 $SUPERSMART_HOME/data/phylota.sqlite 'delete from nodes_194 where common_name like "$SIMDIR%"'
