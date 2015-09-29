#!/bin/bash

# This shell script demonstrates the steps of the SUPERSMART pipeline
# from start to finish, as applied to the phylogeny of the large cats using 
# Lynx as an outgroup. Below we define the names of in- and outgroup and
# a file with a calibration point for the age of the family Felidae

SUBFAMILY=Pantherinae
OUTGROUP=Lynx
FOSSILS=fossils.tsv

# Perform taxonomic name reconciliation on an input list of names.
# Creates a table of NCBI taxonomy identifiers (the taxa table).
# The flag --binomials_only filters out any taxon names that
# are not bonomials, e.g. 'Felis sp. NG192'. The command will by
# default produce a file 'species.tsv' which contains the taxa table.
 
smrt taxize --root_taxa $SUBFAMILY,$OUTGROUP --binomials_only 

# align all phylota clusters for the species in the taxa table.
# produces many aligned fasta files and a file listing these
smrt align --infile species.tsv

# assign orthology among the aligned clusters by reciprocal BLAST
smrt orthologize --infile aligned.txt

# merge the orthologous clusters into a supermatrix with exemplar
# species, two per genus
smrt bbmerge --alnfile merged.txt --taxafile species.tsv

# run an exabayes search on the supermatrix, resulting in a backbone
# posterior sample
smrt bbinfer --supermatrix supermatrix.phy --inferencetool exabayes --cleanup


# root the backbone sample on the outgroup
smrt bbreroot -g $OUTGROUP --backbone backbone.dnd --taxafile species.tsv --smooth

# calibrate the re-rooted backbone tree using treePL
smrt bbcalibrate --tree backbone-rerooted.dnd --supermatrix supermatrix.phy -f $FOSSILS

# build a consensus
smrt consense --infile chronogram.dnd --prob 

# decompose the backbone tree into monophyletic clades. writes a directory
# with suitable alignments for each clade
smrt bbdecompose --backbone consensus.nex --alnfile aligned.txt --taxafile species.tsv

# merge all the alignments for each clades into a nexml file
smrt clademerge --enrich

# run *BEAST for each clade
smrt cladeinfer --ngens=15_000_000 --sfreq=1000 --lfreq=1000

# graft the *BEAST results on the backbone
smrt cladegraft
