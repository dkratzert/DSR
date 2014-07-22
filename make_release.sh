#/bin/bash

GIT="/home/daniel/Downloads/DSR"
PACK="/usr/src/packages"
cd $GIT

dos2unix *

VERSION=$(cat $GIT/dsr.py|grep -e "VERSION ="|cut -d ' ' -f3|tr -d "\'")


cd $GIT

cd setup
sh build_linux_distrib.sh

cd $GIT

cp $GIT/setup/dsr-linux.spec $PACK/SPECS

cd $PACK/SPECS

rpmbuild -bb dsr-linux.spec

cd $GIT

cp -r setup/debian-package/DEBIAN $PACK/BUILD/dsr

cd $PACK/BUILD

dpkg-deb --build dsr

mv dsr.deb dsr-$VERSION.deb 

cp /usr/src/packages/RPMS/noarch/DSR-$VERSION-0.noarch.rpm /usr/src/packages/BUILD/
cp /usr/src/packages/SOURCES/DSR-$VERSION.tar.gz /usr/src/packages/BUILD/
