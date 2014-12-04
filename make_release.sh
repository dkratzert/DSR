#/bin/bash

GIT="/home/daniel/Downloads/DSR"
PACK="/usr/src/packages"
cd $GIT

#dos2unix *

VERSION=$(cat $GIT/dsr.py|grep -e "VERSION ="|cut -d ' ' -f3|tr -d "\'")


cd $GIT

cd setup
sh build_linux_distrib.sh

cd $GIT

cp $GIT/setup/dsr-linux.spec $PACK/SPECS

cd $PACK/SPECS

dos2unix dsr-linux.spec

echo -e "\n####### building rpm #########" 

#rpmbuild -bp dsr-linux.spec -v
rpmbuild -bb dsr-linux.spec -v

echo -e "####### rpm ready #########\n" 

cp /usr/src/packages/RPMS/noarch/DSR-$VERSION-0.noarch.rpm /usr/src/packages/BUILD/
cp /usr/src/packages/SOURCES/DSR-$VERSION.tar.gz /usr/src/packages/BUILD/

#exit 
#######################################
## Debian:
##########################################

cd $GIT

cp -r setup/debian-package/DEBIAN $PACK/BUILD/dsr

cd $PACK/BUILD

dpkg-deb --build dsr

mv dsr.deb dsr-$VERSION.deb 


