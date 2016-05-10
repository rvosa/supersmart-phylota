# This script will automatically generate a debian package from the 
# current SUPERSMART source

# get version
DEBVERSION=$(perl -I$SUPERSMART_HOME/lib -MBio::SUPERSMART -e '($v=$Bio::SUPERSMART::VERSION)=~s/v//; print $v;')

# set up directories
DEBFOLDER=supersmart
DEBFOLDERNAME=$DEBFOLDER-$DEBVERSION
UPSTREAM_TARBALL=supersmart_$DEBVERSION.orig.tar.gz

# make upstream tarball
cp -r upstream/ $DEBFOLDERNAME
tar cvfz $UPSTREAM_TARBALL $DEBFOLDERNAME

# unpack upstream tarball
tar xf $UPSTREAM_TARBALL

cd $DEBFOLDERNAME

# Create the packaging skeleton (debian/*)
dh_make -s -e hannes.hettling@gmail.com --indep

# Remove make calls
#grep -v makefile debian/rules > debian/rules.new 
#mv debian/rules.new debian/rules 

# add docker dependency
sed -i 's/Depends: /Depends: docker-engine (>= 1.9.1-0~), /' debian/control

# add project website
printf "Homepage: http://www.supersmart-project.org\n" >> debian/control

# debian/install must contain the list of scripts to install 
# as well as the target directory
for i in $(ls smrt*); do echo $i usr/bin >> debian/install; done

# Remove the example files
rm debian/*.ex

# make post-installation file which downloads the docker images
printf %"s\n" \
	'#!bin/sh' \
	'echo "Retrieving SUPERSMART image from Docker hub, this might take a while."' \
	'docker pull naturalis/supersmart' > debian/postinst

# Build the package.
# You  will get a lot of warnings and ../somescripts_0.1-1_i386.deb
debuild -us -uc
