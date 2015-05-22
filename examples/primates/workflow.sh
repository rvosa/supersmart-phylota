#!/bin/bash

# this shell script demonstrates the steps of the SUPERSMART pipeline
# from start to finish, as applied to the phylogeny of the primates.

INPUT_NAMES=names.txt
INPUT_FOSSILS=fossils.tsv
OUTGROUP_TAXA=Strepsirrhini

# perform taxonomic name reconciliation on an input list of names.
# creates a table of NCBI taxonomy identifiers (the taxa table).
smrt taxize -i $INPUT_NAMES

# extract the common classification tree for the species in the
# taxa table. produces a newick file with unbranched internal nodes.
smrt classify

# align all phylota clusters for the species in the taxa table.
# produces many aligned fasta files and a file listing these
smrt align 

# assign orthology among the aligned clusters by reciprocal BLAST
smrt orthologize

# merge the orthologous clusters into a supermatrix with exemplar
# species, two per genus
smrt bbmerge

# run an examl search on the supermatrix, resulting in a backbone tree
smrt bbinfer -b 100 -i examl

# root the backbone tree by minimizing the number of non-monophyletic
# higher taxa
smrt bbreroot -g $OUTGROUP_TAXA

# calibrate the re-rooted backbone tree using treePL
smrt bbcalibrate -f $INPUT_FOSSILS

# build a consensus
smrt consense -b 0.0

# decompose the backbone tree into monophyletic clades. writes a directory
# with suitable alignments for each clade
smrt bbdecompose

# merge all the alignments for each clades into a nexml file
smrt clademerge

# run *BEAST for each clade
smrt cladeinfer

# graft the *BEAST results on the backbone
smrt cladegraft
