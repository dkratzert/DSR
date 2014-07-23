#!/bin/bash


TMP="/tmp"
#GIT="/cygdrive/c/Users/$USERNAME/Documents/GitHub/DSR"
#GIT="/home/daniel/GitHub/DSR"
GIT="/home/daniel/Downloads/DSR"
VERSION=$(cat $GIT/dsr.py|grep -e "VERSION ="|cut -d ' ' -f3|tr -d "\'")
BUILDDIR="/usr/src/packages/BUILD/DSR-$VERSION"

echo $VERSION
OUTFILE=DSR-$VERSION.tar

FILES="afix.py
atomhandling.py
atoms.py
constants.py
dbfile.py
dsr.py
dsrparse.py
export.py
misc.py
options.py
resfile.py
restraints.py
terminalsize.py
resi.py
refine.py
pyperclip.py
dsr_db.txt
manuals/DSR-manual.pdf
setup/OlexDSR.py
setup/custom.xld
setup/dsr.sh
dsr
example/p21c.hkl
example/p21c.res
example/p21c_step0.res
example/p21c_step1.res
example/p21c_step2.res
example/p21c_step3.res
example/p21c-step2.ins
"


echo "erstelle tarfile"
rm $GIT/setup/Output/$OUTFILE.gz
touch $GIT/setup/Output/$OUTFILE
#tar cf $GIT/setup/Output/$OUTFILE

TMPDIR=$TMP/DSR-$VERSION
rm -r $TMPDIR
mkdir $TMPDIR
mkdir $TMPDIR/example
mkdir $TMPDIR/setup
mkdir $TMPDIR/manuals

for i in $FILES
do
#    cd $GIT/setup
    cp ../$i $TMPDIR/$i
    dos2unix $TMPDIR/$i -q
    echo "packe " $i
done
# copy networkx
cp -r ../networkx $TMPDIR/

cd $TMP
tar -rf $GIT/setup/Output/$OUTFILE DSR-$VERSION 2> /dev/null

gzip -f $GIT/setup/Output/$OUTFILE

rm -r $BUILDDIR
cp -r $TMPDIR/* $BUILDDIR
# for debian:
DEBDIR="/usr/src/packages/BUILD/dsr"
mkdir $DEBDIR/etc
mkdir $DEBDIR/etc/profile.d
mkdir $DEBDIR/opt
mkdir $DEBDIR/opt/DSR
cp $TMPDIR/* $DEBDIR/opt/DSR
cp -r $TMPDIR/example $DEBDIR/opt/DSR
cp -r $TMPDIR/manuals $DEBDIR/opt/DSR
cp -r $TMPDIR/networkx $DEBDIR/opt/DSR
cp $TMPDIR/setup/* $DEBDIR/etc/profile.d



if [ -e $GIT/setup/Output/$OUTFILE.gz ]
then
    echo "fertig mit Version $VERSION"
    cp $GIT/setup/Output/$OUTFILE.gz /usr/src/packages/SOURCES/
else
    sleep 5
fi
