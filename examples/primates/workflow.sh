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

# run an exabayes search on the supermatrix, resulting in a backbone 
# posterior sample
export SUPERSMART_EXABAYES_NUMGENS="100000"
smrt bbinfer --inferencetool=exabayes --cleanup

# root the backbone sample  on the outgroup
smrt bbreroot -g $OUTGROUP_TAXA

# calibrate the re-rooted backbone tree using treePL
smrt bbcalibrate -f $INPUT_FOSSILS

# clean up temp files from treePL
rm -rf /tmp/*

# build a consensus
smrt consense -b 0.2

# decompose the backbone tree into monophyletic clades. writes a directory
# with suitable alignments for each clade
smrt bbdecompose

# merge all the alignments for each clades into a nexml file
smrt clademerge

# run *BEAST for each clade
smrt cladeinfer --ngens=10_000_000

# graft the *BEAST results on the backbone
smrt cladegraft
