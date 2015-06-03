#!/bin/bash

# This script contains all the steps to create a species tree for European butterflies,
# as selected by the Dutch butterfly foundation ("Vlinderstiching").

# We should have two input files: a list of names, and a fossil table. So far we don't
# have good calibration points but hopefully that will change.
NAMES=names.txt
FOSSILS=fossils.tsv

# Step 1: match the names to the NCBI taxonomy. This is needed because ultimately all
# sequences that the pipeline uses are annotated with NCBI taxonomy identifiers.
smrt taxize --infile=$NAMES

# Step 2: write PhyLoTA alignments that cover the species we have matched in step 1. This
# creates numerous aligned FASTA files and a text file ('aligned.txt') that contains a
# listing of all the alignments. The FASTA files that are created in this step have names
# that are composed as <seed gi>-<seed taxon>-<cluster id>-<type>.fa, meaning:
# - seed gi:    the identifier of the sequence around which this PhyLoTA cluster was built
# - seed taxon: the identifier of the taxon to which <seed gi> belongs
# - cluster id: the internal identifier assigned by PhyLoTA to this cluster
# - type:       the type of cluster, usually a 'subtree' from the NCBI taxonomy
smrt align

# Step 3: write the NCBI taxonomy as a classification tree. This creates a newick tree
# file ('classification-tree.dnd') with additional, unbranched interior nodes for all
# applicable taxonomic levels. The tree will likely be highly polytomous, so in steps
# where we use it as a starting tree we randomly resolve it, and collapse all the
# unbranched interior nodes.
smrt classify

# Step 4: assess orthology among PhyLoTA alignments. This step is needed because PhyLoTA
# is unable to perform all-vs-all BLASTing across the entire GenBank (consider the
# combinatorics here) so instead we do this on the fly just for the alignments that were
# created in step 2. This step will create numerous aligned FASTA files, each called
# clusterXXX.fa, where XXX is an integer. In addition, a simple text file ('merged.txt')
# will contain a listing of these alignments.
smrt orthologize

# Step 5: build a supermatrix. This is done by concatenating selected orthologous clusters
# from Step 4. The clusters are selected by a number of interacting criteria, which can
# be altered in the configuration file 'supersmart.ini'
# - BACKBONE_MAX_DISTANCE: the maximum average pairwise distance below which we still
#                          accept the cluster. Above this value we assume that the cluster
#                          is too saturated to be of use.
# - BACKBONE_MIN_COVERAGE: the minimum number of times a taxon must participate in loci
#                          in the supermatrix in order to include it.
# We can override these config values from the environment, like so:
export SUPERSMART_BACKBONE_MAX_DISTANCE="0.05"
export SUPERSMART_BACKBONE_MIN_COVERAGE="3"
smrt bbmerge

# Step 6: build a backbone tree. There are multiple inference tools that can be used for
# this (raxml, exabayes, examl). Here we use exabayes. The inference is done on the
# supermatrix constructed in step 5, using the classification tree from step 3 as a
# starting tree. In this example we follow EXABAYES defaults from supersmart.ini. As a
# result, the output file ('backbone.dnd') will be a newick tree file a sample from the
# posterior distribution in it. This step normally creates very many intermediate files
# (bootstrapped matrices, various log and checkpoint files produced by the inference tool,
# separate output tree files), which we clean up by providing the '--cleanup' flag.
export SUPERSMART_EXABAYES_NUMGENS="100000"
smrt bbinfer --inferencetool=exabayes --cleanup

# Step 7: reroot the backbone trees. The trees resulting from step 6 are unrooted. There
# are different ways to root these (using outgroups or by picking the rooting that best
# fits the taxonomy). Here we fit to the taxonomy: we have two genera so this amounts to
# the same thing as considering either of these the outgroup with respect to the other.
# By default produces a file 'backbone-rerooted.dnd'
smrt bbreroot --smooth

# Step 8: calibrate the backbone trees from step 7. This step uses treePL to create
# ultrametric trees ('chronogram.dnd') using a penalized likelihood approach with
# age constraints as extracted from the fossils file.
smrt bbcalibrate --fossiltable=$FOSSILS

# Step 9: build a consensus tree. As we ran exabayes in step 6 we will want to discard
# a burnin. This would be different had we done a bootstrapping analysis in step 6.
smrt consense --burnin=0.20 --prob

# Step 10: decompose the backbone into clades. This step traverses the consensus tree from
# step 9 and breaks it up into monophyletic clades (in principle, genera, unless these are
# not monophyletic, in which case the decomposition happens at the nearest ancestor that
# subtends a monophyletic group of genera). For each clade, alignments (produced in step
# 2) are selected that meet the following criteria:
# - CLADE_MAX_DISTANCE: maximum average pairwise distance to accept the alignment
# - CLADE_MIN_DENSITY:  minimum number of clade members that must be present
export SUPERSMART_CLADE_MAX_DISTANCE="0.1"
export SUPERSMART_CLADE_MIN_DENSITY="0.3"
smrt bbdecompose

# Step 11: for each clade, merge the separate clade alignments from step 10 into a single
# input NeXML file for *BEAST.
smrt clademerge --enrich

# Step 12: for each clade, run *BEAST. By default this uses a very small number of
# generations (100,000), which is strictly intended for testing. It is not possible to
# make any general recommendations for what the right number is because this depends on a
# lot of different variables (number of taxa, signal in the data, mixing of the chains)
# so for publishable results convergence should be checked, for example using 'tracer'.
# For this specific example, 20,000,000 generations for each clade appears to work.
smrt cladeinfer --ngens=20_000_000 --sfreq=1000 --lfreq=1000

# Step 13: graft the clade trees onto the backbone.
# We need some general solutions for negative branch lengths. These can happen for
# example when a pair of exemplars actually doesn't cross the root of the clade they
# represent. Scaling on their MRCA very often results in the actual root of the clade
# being older than its parent, which in turn results in negative branch lengths. Some
# approaches:
# 1. set negative branch lengths to zero
# 2. rescale the whole clade (+ its root branch) to fit below its parent
# 3. push the parent deeper (might just move the problem)
smrt cladegraft

# Step 14: plot the final result
smrt-utils plot
