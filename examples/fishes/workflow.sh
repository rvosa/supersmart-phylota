app/smrt taxize -w examples/fishes/ -i examples/fishes/names.txt

app/smrt classify -w examples/fishes/ -i examples/fishes/species.tsv

app/smrt align -w examples/fishes/ -i examples/fishes/species.tsv

app/smrt orthologize -w examples/fishes/ -i examples/fishes/aligned.txt

app/smrt bbmerge -w examples/fishes/ -a examples/fishes/merged.txt -t examples/fishes/species.tsv

app/smrt bbinfer -w examples/fishes/ -s examples/fishes/supermatrix.phy -t examples/fishes/classification-tree.dnd


