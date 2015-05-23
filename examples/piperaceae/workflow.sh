#!/bin/bash

# This script demonstrate all the steps to create a tree for species from the genera Piper 
# and Peperomia that occur in the Guyanas, i.e. a small case where we have a set of names 
# that were selected somehow (in this case, by text mining the Flora of the Guyanas) and 
# for which we want to infer a time-calibrated tree.

# We have two input files: a list of names, and a fossil table. The fossil table contains
# a single calibration point based on a statement in the literature that the genus Piper's
# crown age is ~30MY. For all the intermediate files that will be created we will use the
# default names that the pipeline expects to be present in the working directory. 
# Alternative names and locations can of course be provided to each command; consult the
# help messages of each command to learn how to do this. With the settings defined here
# (i.e. local names inside a folder) the subsequent commands need to be executed within
# this same folder.
export NAMES=names.txt
export FOSSILS=fossils.tsv
export OUTGROUP=Peperomia

# Step 1: match the names to the NCBI taxonomy. This is needed because ultimately all 
# sequences that the pipeline uses are annotated with NCBI taxonomy identifiers.
smrt taxize -i $NAMES

# Step 2: write PhyLoTA alignments that cover the species we have matched in step 1. This
# creates numerous aligned FASTA files and a text file ('aligned.txt') that contains a
# listing of all the alignments. The FASTA files that are created in this step have names
# that are composed as <seed gi>-<seed taxon>-<cluster id>-<type>.fa, meaning:
# - seed gi:    the identifier of the sequence around which this PhyLoTA cluster was built
# - seed taxon: the identifier of the taxon to which <seed gi> belongs
# - cluster id: the internal identifier assigned by PhyLoTA to this cluster
# - type:       the type of cluster, usually a 'subtree' from the NCBI taxonomy
smrt align -d alignments

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
smrt bbmerge

# Step 6: build a backbone tree. There are multiple inference tools that can be used for
# this (raxml, exabayes, examl). The default is examl. The inference is done on the 
# supermatrix constructed in step 5, using the classification tree from step 3 as a 
# starting tree. In this example we will run 100 bootstrap replicates. As a result, the
# output file ('backbone.dnd') will be a newick tree file with 100 ML trees in it. This
# step normally creates very many intermediate files (bootstrapped matrices, various
# log and checkpoint files produced by the inference tool, separate output tree files),
# which we clean up by providing the '-x' flag.
smrt bbinfer -b 100 -x

# Step 7: reroot the backbone trees. The trees resulting from step 6 are unrooted. There
# are different ways to root these (using outgroups or by picking the rooting that best 
# fits the taxonomy). We will root on the group Piper
smrt bbreroot -g $OUTGROUP

# Step 8: calibrate the backbone trees from step 7. This step uses treePL to create 
# ultrametric trees ('chronogram.dnd') using a penalized likelihood approach with 
# age constraints as extracted from the fossils file.
smrt bbcalibrate -f $FOSSILS

# Step 9: build a consensus tree. In this example we have 100 bootstrap trees that we have
# rerooted and calibrated. By now building a consensus we will end up with a single tree
# that will have bootstrap values and bootstrapped age ranges on the branches. The 
# consensus is calculated using treeannotator, and so the metadata that is attached to the
# nodes can be viewed by opening the resulting file ('consensus.nex') in figtree. As the 
# trees are bootstrap trees we don't throw out a burnin fraction. This would be different 
# had we used exabayes in step 6.
smrt consense --burnin=0

# Step 10: decompose the backbone into clades. This step traverses the consensus tree from
# step 9 and breaks it up into monophyletic clades (in principle, genera, unless these are
# not monophyletic, in which case the decomposition happens at the nearest ancestor that
# subtends a monophyletic group of genera). For each clade, alignments (produced in step 
# 2) are selected that meet the following criteria:
# - CLADE_MAX_DISTANCE: maximum average pairwise distance to accept the alignment
# - CLADE_MIN_DENSITY:  minimum number of clade members that must be present
smrt bbdecompose

# Step 11: for each clade, merge the separate clade alignments from step 10 into a single
# input NeXML file for *BEAST.
smrt clademerge

# Step 12: for each clade, run *BEAST. By default this uses a very small number of 
# generations (100,000), which is strictly intended for testing. It is not possible to 
# make any general recommendations for what the right number is because this depends on a
# lot of different variables (number of taxa, signal in the data, mixing of the chains)
# so for publishable results convergence should be checked, for example using 'tracer'.
# For this specific example, 30,000,000 generations for each clade appears to work.
smrt cladeinfer -n 30_000_000

# Step 13: graft the clade trees onto the backbone.
# XXX ANNOTATIONS LOST WHEN READING IN CONSENSE TREE!
smrt cladegraft

# Step 14: plot the final result
smrt-utils plot
