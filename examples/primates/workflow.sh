smrt taxize -i names.txt
smrt classify -i species.tsv
smrt align -i species.tsv
smrt orthologize -i aligned.txt
smrt bbmerge -a merged.txt-t species.tsv
smrt bbinfer -s supermatrix.phy-t classification-tree.dnd-i examl
smrt bbreroot -b backbone.dnd-t species.tsv
smrt bbcalibrate -t backbone-rerooted.dnd-s supermatrix.phy-f fossils.tsv-o backbone-calibrated.dnd
smrt bbdecompose -b backbone-calibrated.dnd-c classification-tree.dnd-a merged.txt-t species.tsv
smrt clademerge  
smrt cladeinfer 
smrt cladegraft -b backbone-calibrated.dnd
