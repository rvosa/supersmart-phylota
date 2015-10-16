#!/bin/bash

# This script makes synthetic datasets and runs the
# SUPERSMART pipline on them. Prerequisite is a 
# finalized run of the SUPERSMART pipeline.

# number of replicated
REPLICATES="3"

# prefix for simulation directories
PREFIX="sim-"

# taxon set parameters
export OUTGROUP="Lynx"

# backbone parameters
export SUPERSMART_BACKBONE_MAX_DISTANCE="0.2"
export SUPERSMART_BACKBONE_MIN_COVERAGE="3"
export SUPERSMART_BACKBONE_MAX_COVERAGE="5"
export SUPERSMART_EXABAYES_NUMGENS="100000"
export SUPERSMART_MSA_TOOL="mafft"

# clade parameters
export SUPERSMART_CLADE_MAX_DISTANCE="0.1"
export SUPERSMART_CLADE_MIN_DENSITY="0.5"
export SUPERSMART_CLADE_TAXON_MIN_MARKERS="2"
export SUPERSMART_CLADE_MAX_MARKERS="10"

# make synthetic datasets
for i in $( seq 1 $REPLICATES ); do
	sh make-replicated-dataset.sh $PREFIX$i
done;

# make plots of original and synthetic datasets
if [ ! -e "alignment-stats.tsv" ]; then
	smrt-utils alnstats -a aligned.txt
fi
ARGS="alignment-stats.tsv"
for file in $( ls sim*/alignment-stats.tsv ); do
	ARGS="$ARGS $file"
done
Rscript plot-aln-stats.R $ARGS

# run SUPERSMART on synthetic datasets 
for i in $( seq 1 $REPLICATES ); do
	sh process-replicated-dataset.sh $PREFIX$i
done;
