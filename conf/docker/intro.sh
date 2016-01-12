#!/bin/bash

# make custom prompt with version and commit id 
cd $SUPERSMART_HOME
commit=`git log -1 | head -1 | cut -f 2 -d ' '`
hash=${commit:0:10}
cd - >/dev/null
VERSION=`perl -MBio::SUPERSMART -e 'print $Bio::SUPERSMART::VERSION'`
export PS1="\[\e[0;32m\]SUPERSMART $VERSION - $hash\[\e[m\] \[\e[1;34m\]\w\[\e[m\] \[\e[1;32m\]\$\[\e[m\] "

# unzip phylota on first use
if [ -e $SUPERSMART_HOME/data/phylota.sqlite.gz ]; then
	echo 'preparing database for usage'
	gunzip $SUPERSMART_HOME/data/phylota.sqlite.gz
fi

# print welcome message
echo '
   _____ _    _ _____  ______ _____    
  / ____| |  | |  __ \|  ____|  __ \   
 | (___ | |  | | |__) | |__  | |__) |  
  \___ \| |  | |  ___/|  __| |  _  /   
  ____) | |__| | |    | |____| | \ \   
 |_____/ \____/|_|    |______|_|__\_\_ 
  / ____|  \/  |   /\   |  __ \__   __|
 | (___ | \  / |  /  \  | |__) | | |   
  \___ \| |\/| | / /\ \ |  _  /  | |   
  ____) | |  | |/ ____ \| | \ \  | |   
 |_____/|_|  |_/_/    \_\_|  \_\ |_|   
'
echo "Welcome to the SUPERSMART pipeline developed by 
Naturalis Biodiversity Centre & Gothenburg University.

Available commands:

     smrt
     smrt-utils
     smrt-config 

To get started, try 'smrt --help' and visit www.supersmart-project.org.
"
