#!/bin/bash

# this shell script demonstrates the steps of the SUPERSMART pipeline
# from start to finish, as applied to the phylogeny of the primates.

# perform taxonomic name reconciliation on an input list of names.
# creates a table of NCBI taxonomy identifiers (the taxa table).
smrt taxize -i names.txt

# extract the common classification tree for the species in the
# taxa table. produces a newick file with unbranched internal nodes.
smrt classify -i species.tsv

# align all phylota clusters for the species in the taxa table.
# produces many aligned fasta files and a file listing these
smrt align -i species.tsv

# assign orthology among the aligned clusters by reciprocal BLAST
smrt orthologize -i aligned.txt

# merge the orthologous clusters into a supermatrix with exemplar
# species, two per genus
smrt bbmerge -a merged.txt -t species.tsv

# run an examl search on the supermatrix, resulting in a backbone tree
smrt bbinfer -s supermatrix.phy -t classification-tree.dnd -i examl

# root the backbone tree by minimizing the number of non-monophyletic
# higher taxa
smrt bbreroot -b backbone.dnd -t species.tsv

# calibrate the re-rooted backbone tree using treePL
smrt bbcalibrate -t backbone-rerooted.dnd -s supermatrix.phy -f fossils.tsv -o backbone-calibrated.dnd

# decompose the backbone tree into monophyletic clades. writes a directory
# with suitable alignments for each clade
smrt bbdecompose -b backbone-calibrated.dnd -c classification-tree.dnd -a aligned.txt -t species.tsv

# merge all the alignments for each clades into a nexml file
smrt clademerge

# run *BEAST for each clade
smrt cladeinfer

# graft the *BEAST results on the backbone
smrt cladegraft -b backbone-calibrated.dnd
