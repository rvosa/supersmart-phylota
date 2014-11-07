app/smrt taxize -w examples/primates/ -i examples/primates/names.txt

app/smrt classify -w examples/primates/ -i examples/primates/species.tsv

app/smrt align -w examples/primates/ -i examples/primates/species.tsv

app/smrt orthologize -w examples/primates/ -i examples/primates/aligned.txt

app/smrt bbmerge -w examples/primates/ -a examples/primates/merged.txt -t examples/primates/species.tsv

app/smrt bbinfer -w examples/primates/ -s examples/primates/supermatrix.phy -t examples/primates/classification-tree.dnd -i examl

app/smrt bbreroot -w examples/primates/ -b examples/primates/backbone.dnd -t examples/primates/species.tsv

app/smrt bbcalibrate -t examples/primates/backbone-rerooted.dnd -s examples/primates/supermatrix.phy -f examples/primates/fossils.tsv -o examples/primates/backbone-calibrated.dnd

app/smrt bbdecompose -w examples/primates/ -b examples/primates/backbone-calibrated.dnd -c examples/primates/classification-tree.dnd -a examples/primates/merged.txt -t examples/primates/species.tsv

app/smrt clademerge -w examples/primates/ 

app/smrt cladeinfer -w examples/primates

app/smrt cladegraft -w examples/primates/ -b examples/primates/backbone-calibrated.dnd 
