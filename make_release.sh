#/bin/bash

GIT="/home/daniel/Downloads/DSR"
PACK="/usr/src/packages"

VERSION=$(cat $GIT/dsr.py|grep -e "VERSION ="|cut -d ' ' -f3|tr -d "\'")


cd $GIT

git pull

cp $GIT/setup/dsr-linux.spec $PACK/SPECS

cd $PACK/SPECS

rpmbuild -bb dsr-linux.spec

sh setup/build_linux_distrib.sh

cp -r setup/debian-package/DEBIAN $PACK/BUILD/dsr

cd $PACK/BUILD

dpkg-deb --build dsr

mv dsr.deb dsr-$VERSION.deb 

