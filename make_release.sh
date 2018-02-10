#!/usr/bin/env bash
#/bin/bash

GIT="$(pwd)"
PACK="/usr/src/packages"

cd $GIT

VERSION=$(cat ${GIT}/dsr.py|grep -e "VERSION ="|cut -d ' ' -f3|tr -d "\'")
TMPDIR=/tmp/DSR-${VERSION}
BUILDDIR=/usr/src/packages/BUILD/DSR-${VERSION}
DEBDIR=/usr/src/packages/BUILD/dsr

mkdir ${BUILDDIR}

cd $GIT

#cd setup
#sh ./setup/build_linux_distrib.sh
rm $GIT/setup/Output/DSR-${VERSION}.tar.gz
PYTHONPATH=$GIT python ./scripts/make_zipfile.py
cp setup/Output/DSR-${VERSION}.tar.gz /usr/src/packages/BUILD/

cd ${GIT}

############  Prepare debian package  ########################

rm -rf /tmp/DSR-${VERSION}
tar -xf setup/Output/DSR-$VERSION.tar.gz -C /tmp
cp -r $TMPDIR/* $BUILDDIR
# for debian:

rm -r ${DEBDIR}
mkdir $DEBDIR
mkdir $DEBDIR/etc
mkdir $DEBDIR/etc/profile.d
mkdir $DEBDIR/opt
mkdir $DEBDIR/opt/DSR
cp $TMPDIR/* $DEBDIR/opt/DSR
cp -r $TMPDIR/example $DEBDIR/opt/DSR
cp -r $TMPDIR/manuals $DEBDIR/opt/DSR
cp -r $TMPDIR/networkx $DEBDIR/opt/DSR
cp -r $TMPDIR/mpmath $DEBDIR/opt/DSR
cp $TMPDIR/setup/* $DEBDIR/etc/profile.d

##########################################

cp $GIT/setup/dsr-linux.spec $PACK/SPECS

cd $PACK/SPECS

dos2unix dsr-linux.spec

echo -e "\n####### building rpm #########" 

#rpmbuild -bp dsr-linux.spec -v
rpmbuild -bb dsr-linux.spec -v

echo -e "####### rpm ready #########\n" 

cp /usr/src/packages/RPMS/noarch/DSR-$VERSION-0.noarch.rpm /usr/src/packages/BUILD/

#exit 
#######################################
## Debian:
##########################################

cd $GIT

cp -r setup/debian-package/DEBIAN $PACK/BUILD/dsr

cd $PACK/BUILD

dpkg-deb --build dsr

mv dsr.deb dsr-$VERSION.deb 


