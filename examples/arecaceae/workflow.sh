#!/bin/bash

FOSSILS=fossils.tsv
OUTGROUP=Anigozanthos,Catopsis,Typha
SPECIES=species.tsv

# clear database
 sqlite3 $SUPERSMART_HOME/data/phylota.sqlite 'delete from seqs where acc like "SMRT%";'
 sqlite3 $SUPERSMART_HOME/data/phylota.sqlite 'delete from nodes_194 where common_name like "SMRT%";'

# insert species that are not yet in database
 smrt-utils dbinsert -t $SPECIES

# insert custom sequences into database
# this produces output files with suffix '-smrt-inserted.fa' which
# will then be appended to the alignment list
 smrt-utils dbinsert -a matK.phy -f phylip
 smrt-utils dbinsert -a rpb2.phy -f phylip


# align all phylota clusters for the species in the taxa table.
# produces many aligned fasta files and a file listing these
export SUPERSMART_MSA_TOOL="mafft"
smrt align

for i in $(ls *-smrt-inserted.fa);do echo $PWD/$i; done >> aligned.txt

# assign orthology among the aligned clusters by reciprocal BLAST
export SUPERSMART_BACKBONE_MAX_DISTANCE="0.20"
smrt orthologize

# merge the orthologous clusters into a supermatrix with exemplar
# species, two per genus
export SUPERSMART_BACKBONE_MIN_COVERAGE="3"
export SUPERSMART_BACKBONE_MAX_COVERAGE="5"
smrt bbmerge -o supermatrix.phy

# run an exabayes search on the supermatrix, resulting in a backbone
# posterior sample
export SUPERSMART_EXABAYES_NUMGENS="100000"
export SUPERSMART_EXABAYES_NUMRUNS="4"
export SUPERSMART_EXABAYES_NUMCHAINS="2"
smrt bbinfer --inferencetool=exabayes -s supermatrix.phy -o backbone-2.dnd --cleanup

# root the backbone sample  on the outgroup
smrt bbreroot -g $OUTGROUP --smooth

# calibrate the re-rooted backbone tree using treePL
smrt bbcalibrate -f $FOSSILS

# build a consensus
smrt consense -b 0.1 --prob

# decompose the backbone tree into monophyletic clades. writes a directory
# with suitable alignments for each clade
export SUPERSMART_CLADE_MAX_DISTANCE="0.1"
export SUPERSMART_CLADE_MIN_DENSITY="0.3"
export SUPERSMART_CLADE_MIN_COVERAGE="2"
export SUPERSMART_CLADE_MAX_COVERAGE="10" 
smrt bbdecompose

# merge all the alignments for each clades into a nexml file
smrt clademerge --enrich

# run *BEAST for each clade
smrt cladeinfer --ngens=20_000_000 --sfreq=1000 --lfreq=1000

# graft the *BEAST results on the backbone
smrt cladegraft
