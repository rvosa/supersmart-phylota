## INITIAL TREE INFERENCE

FOSSILS=fossils.tsv
smrt taxize -r Coturnix,Odontophoridae
smrt align
smrt classify
smrt orthologize
export SUPERSMART_BACKBONE_MAX_DISTANCE=0.2
export SUPERSMART_BACKBONE_MIN_COVERAGE=3
smrt bbmerge -v
smrt bbinfer --cleanup
smrt bbreroot --smooth
smrt bbcalibrate --fossiltable=$FOSSILS
smrt consense --burnin=0.20 --prob
export SUPERSMART_CLADE_MAX_DISTANCE=0.3
export SUPERSMART_CLADE_MIN_DENSITY=0.2
export SUPERSMART_CLADE_MAX_MARKERS=10
export SUPERSMART_CLADE_TAXON_MIN_MARKERS=2
export SUPERSMART_CLADE_MAX_HAPLOTYPES=3
smrt bbdecompose
smrt clademerge 
smrt cladeinfer --ngens=100_000 --sfreq=100 --lfreq=100
smrt cladegraft

mkdir simulations
cp final.nex simulations
cp aligned.txt simulations
cp fossils.tsv simulations
cd simulations

## REPLICATION OF DATASET
smrt-utils replicate -t final.nex -f nexus -a aligned.txt -l replicate.log -v
smrt-utils dbinsert -s aligned-replicated.txt -t taxa-replicated.tsv

## RERUN ANALYSIS WITH REPLICATED DATASET
smrt classify -i taxa-replicated.tsv
smrt orthologize -i aligned-smrt-inserted.txt 
smrt bbmerge -t taxa-replicated.tsv -a merged.txt
smrt bbinfer -b 10 -s supermatrix.phy
smrt bbreroot -b backbone.dnd -t taxa-replicated.tsv
smrt bbcalibrate -t backbone-rerooted.dnd -f fossils.tsv
smrt consense -i chronogram.dnd
smrt bbdecompose -b consensus.nex -a aligned-smrt-inserted.txt -t taxa-replicated.tsv
smrt clademerge --enrich
smrt cladeinfer --ngens=100_000 --sfreq=100 --lfreq=100
smrt cladegraft

