#!/bin/bash

DISTANCES="0.05 0.07 0.09 0.11 0.13 0.15 0.17 0.19 0.21 0.23 0.25"
COVERAGES="1 3 5 7 9 11 13 15 17 19"
OUTGROUP=Strepsirrhini
export SUPERSMART_EXABAYES_NUMGENS="100000"

for DIST in $DISTANCES; do
	for COVER in $COVERAGES; do

		# file name of supermatrix
                SUPERMATRIX="supermatrix-$DIST-$COVER.phy"

		# skip over ones we've already done
		if [ -e "$SUPERMATRIX" ]; then
			echo "already done $SUPERMATRIX"
			continue
		fi

		# write current iteration to STDOUT/nohup.out
		echo "dist: $DIST - cover: $COVER"

		# remove previous superclusters, if any
		rm cluster*.fa
		rm seeds.fa*

		# set environment variables
		export SUPERSMART_BACKBONE_MAX_DISTANCE=$DIST
		export SUPERSMART_BACKBONE_MIN_COVERAGE=$COVER

		# make superclusters
		smrt orthologize

		# merge to supermatrix
		smrt bbmerge -o $SUPERMATRIX

		# run the search
		BACKBONE="backbone-$DIST-$COVER.dnd"
		smrt bbinfer --inferencetool=exabayes --cleanup -s $SUPERMATRIX -o $BACKBONE

		# reroot
		BACKBONE_REROOTED="backbone-$DIST-$COVER-rerooted.dnd"
		smrt bbreroot -g $OUTGROUP -b $BACKBONE -o $BACKBONE_REROOTED

		# build consensus
		CONSENSE=consensus-$DIST-$COVER.nex
		smrt consense -b 0.2 --prob -i $BACKBONE_REROOTED -o $CONSENSE

		# store markers table
		mv markers-backbone.tsv markers-backbone-$DIST-$COVER.tsv
	done
done
