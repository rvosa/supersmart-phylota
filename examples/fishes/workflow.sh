app/smrt taxize -w examples/fishes/ -i examples/fishes/names.txt

app/smrt classify -w examples/fishes/ -i examples/fishes/species.tsv

app/smrt align -w examples/fishes/ -i examples/fishes/species.tsv

app/smrt orthologize -w examples/fishes/ -i examples/fishes/aligned.txt

app/smrt bbmerge -w examples/fishes/ -a examples/fishes/merged.txt -t examples/fishes/species.tsv

app/smrt bbinfer -w examples/fishes/ -s examples/fishes/supermatrix.phy -t examples/fishes/classification-tree.dnd -i raxml

app/smrt bbreroot -w examples/fishes/ -b examples/fishes/backbone.dnd -t examples/fishes/species.tsv

app/smrt bbdecompose -w examples/fishes/ -b examples/fishes/backbone-rerooted.dnd -c examples/fishes/classification-tree.dnd -a examples/fishes/merged.txt -t examples/fishes/species.tsv

app/smrt clademerge -w examples/fishes/ 

app/smrt cladeinfer -w examples/fishes

app/smrt cladegraft -w examples/fishes/ -b examples/fishes/backbone-rerooted.dnd 