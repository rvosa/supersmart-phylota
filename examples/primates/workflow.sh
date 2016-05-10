#!/bin/bash

# this shell script demonstrates the steps of the SUPERSMART pipeline
# from start to finish, as applied to the phylogeny of the primates.

FOSSILS=fossils.tsv
OUTGROUP=Strepsirrhini

# perform taxonomic name reconciliation on an input list of names.
# creates a table of NCBI taxonomy identifiers (the taxa table).
smrt taxize -r Primates -b

# align all phylota clusters for the species in the taxa table.
# produces many aligned fasta files and a file listing these
smrt align

# assign orthology among the aligned clusters by reciprocal BLAST
export SUPERSMART_BACKBONE_MAX_DISTANCE="0.1"
smrt orthologize

# merge the orthologous clusters into a supermatrix with exemplar
# species, two per genus
export SUPERSMART_BACKBONE_MIN_COVERAGE="3"
export SUPERSMART_BACKBONE_MAX_COVERAGE="5"
smrt bbmerge

# run an exabayes search on the supermatrix, resulting in a backbone
# posterior sample
export SUPERSMART_EXABAYES_NUMGENS="100000"
smrt bbinfer --inferencetool=exabayes --cleanup

# root the backbone sample  on the outgroup
smrt bbreroot -g $OUTGROUP --smooth

# calibrate the re-rooted backbone tree using treePL
smrt bbcalibrate -f $FOSSILS

# build a consensus
smrt consense -b 0.2 --prob

# decompose the backbone tree into monophyletic clades. writes a directory
# with suitable alignments for each clade
export SUPERSMART_CLADE_MAX_DISTANCE="0.2"
export SUPERSMART_CLADE_MIN_DENSITY="0.5"
export SUPERSMART_CLADE_MIN_COVERAGE="2"
export SUPERSMART_CLADE_MAX_COVERAGE="10"

smrt bbdecompose -g

# merge all the alignments for each clades into a nexml file
smrt clademerge --enrich

# run *BEAST for each clade
smrt cladeinfer --ngens=30_000_000 --sfreq=1000 --lfreq=1000

# graft the *BEAST results on the backbone
smrt cladegraft
