#!/bin/bash

# bash script to update version numbers as per 
# issue #41 before kicking off a packer run

# most recent tag version on git
TAGVERSION=`git describe --tags | cut -f1 -d '-' | head -1`

# version in framework code
APIVERSION=`perl -I../lib -MBio::SUPERSMART -e "print $Bio::SUPERSMART::VERSION"`

if [ "$TAGVERSION" -ne "$APIVERSION" ]; then

	# this so that we see what's happening in the travis logs
	echo "changing from version $APIVERSION to $TAGVERSION"

	# replace in the base module file
	sed -e "s/$APIVERSION/$TAGVERSION/" ../lib/Bio/SUPERSMART.pm > ../lib/Bio/SUPERSMART.pm.bak
	mv ../lib/Bio/SUPERSMART.pm.bak ../lib/Bio/SUPERSMART.pm
	
	# replace in packer template
	APIVERSION=`echo $APIVERSION | sed -e 's/v//'`
	TAGVERSION=`echo $TAGVERSION | sed -e 's/v//'`
	sed -e "s/$APIVERSION/$TAGVERSION/" ../conf/template.json > ../conf/template.json.bak
	mv ../conf/template.json.bak ../conf/template.json
else
	echo "version unchanged: (tag) $TAGVERSION == (api) $APIVERSION"
fi