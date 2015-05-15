#!/bin/bash

# bash script to update version numbers as per 
# issue #41 before kicking off a packer run

# shorter copy of SUPERSMART_HOME
H=$SUPERSMART_HOME

# most recent tag version on git
TAGVERSION=`git describe --tags | cut -f1 -d '-' | head -1`

# version in framework code
APIVERSION=`perl -I$H/lib -MBio::SUPERSMART -e 'print $Bio::SUPERSMART::VERSION'`

# operating system. need to know this because OSX wants an empty string
# for `sed -i`, whereas linux does not
OS=`perl -e 'print $^O'`

if [ "$TAGVERSION" != "$APIVERSION" ]; then

	# this so that we see what's happening in the travis logs
	echo "changing from version $APIVERSION to $TAGVERSION"

	# replace in the base module file
	echo "editing $H/lib/Bio/SUPERSMART.pm in place"
	if [ "$OS" == "darwin" ]; then
		sed -i '' -e "s/$APIVERSION/$TAGVERSION/" $H/lib/Bio/SUPERSMART.pm
	else
		sed -i -e "s/$APIVERSION/$TAGVERSION/" $H/lib/Bio/SUPERSMART.pm	
	fi
	
	# replace in packer template
	echo "editing $H/conf/template.json in place"
	APIVERSION=`echo $APIVERSION | sed -e 's/v//'`
	TAGVERSION=`echo $TAGVERSION | sed -e 's/v//'`
	if [ "$OS" == "darwin" ]; then
		sed -i '' -e "s/$APIVERSION/$TAGVERSION/" $H/conf/template.json
	else
		sed -i -e "s/$APIVERSION/$TAGVERSION/" $H/conf/template.json
	fi
else
	echo "version unchanged: (tag) $TAGVERSION == (api) $APIVERSION"
fi
