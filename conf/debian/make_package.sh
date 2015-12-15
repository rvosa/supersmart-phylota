DEBVERSION=0.1
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

sed -i 's/Depends: /Depends: docker-engine (>= 1.9.1-0~), /' debian/control

# debian/install must contain the list of scripts to install 
# as well as the target directory
for i in $(ls smrt*); do echo $i usr/bin >> debian/install; done

# Remove the example files
rm debian/*.ex

# Build the package.
# You  will get a lot of warnings and ../somescripts_0.1-1_i386.deb
debuild -us -uc
