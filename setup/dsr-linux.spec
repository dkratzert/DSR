Summary: DSR - A program for modelling of disordered solvents with SHELXL
Name: DSR
Provides: DSR
Packager: dkratzert@gmx.de
Version: 1.5.0
Release: 0
Requires: python, xclip
BuildRoot: %{_tmppath}/%{name}-%{version}-build/opt/DSR
BuildArch: noarch
URL: https://www.xs3.uni-freiburg.de/research/dsr
License: Beerware
Source: DSR-%{version}.tar.gz
Group: Productivity/Scientific/Chemistry
%description
This program consists of a text database with fragments of molecules and 
the DSR program. It acts as a preprocessor for SHELXL .res files. The user 
inserts a special command in the SHELXL .res file and the DSR program reads 
this information to put a molecule or fragment with the desired atoms on the 
position of the target atoms specified by the user. Bond restraints are applied
from the database to the molecule.
Development is on GitHub: https://github.com/dkratzert/dsr

%prep
%setup -q
mkdir -p /opt/DSR
mkdir -p /opt/DSR/manuals
mkdir -p /opt/DSR/example
mkdir -p /opt/DSR/setup

%post
if [ ! -f /opt/DSR/dsr_user_db.txt ]
then
    touch /opt/DSR/dsr_user_db.txt
    chmod a+rw /opt/DSR/dsr_user_db.txt
else
    chmod a+rw /opt/DSR/dsr_user_db.txt
fi
chmod a+rw /opt/DSR/example

%install
install -m 644 afix.py /opt/DSR/afix.py
install -m 644 atomhandling.py /opt/DSR/atomhandling.py
install -m 644 atoms.py /opt/DSR/atoms.py
install -m 644 constants.py /opt/DSR/constants.py
install -m 644 dbfile.py /opt/DSR/dbfile.py
install -m 644 dsr.py /opt/DSR/dsr.py
install -m 755 dsr /opt/DSR/dsr
install -m 644 dsrparse.py /opt/DSR/dsrparse.py
install -m 644 export.py /opt/DSR/export.py
install -m 644 misc.py /opt/DSR/misc.py
install -m 644 options.py /opt/DSR/options.py
install -m 644 terminalsize.py /opt/DSR/terminalsize.py
install -m 644 resfile.py /opt/DSR/resfile.py
install -m 644 restraints.py /opt/DSR/restraints.py
install -m 644 resi.py /opt/DSR/resi.py
install -m 644 refine.py /opt/DSR/refine.py
install -m 644 pyperclip.py /opt/DSR/pyperclip.py
install -m 644 dsr_db.txt /opt/DSR/dsr_db.txt
install -m 644 manuals/DSR-manual.pdf /opt/DSR/manuals/DSR-manual.pdf
install -m 644 setup/dsr.sh /etc/profile.d/dsr.sh
install -m 644 setup/OlexDSR.py /opt/DSR/setup/OlexDSR.py
install -m 644 setup/custom.xld /opt/DSR/setup/custom.xld
install -m 644 example/p21c.hkl /opt/DSR/example/p21c.hkl
install -m 644 example/p21c.res /opt/DSR/example/p21c.res
install -m 644 example/p21c_step0.res /opt/DSR/example/p21c_step0.res
install -m 644 example/p21c_step1.res /opt/DSR/example/p21c_step1.res
install -m 644 example/p21c_step2.res /opt/DSR/example/p21c_step2.res
install -m 644 example/p21c_step3.res /opt/DSR/example/p21c_step3.res
install -m 644 example/p21c-step2.ins /opt/DSR/example/p21c-step2.ins
install -d -m 755 networkx /opt/DSR/networkx
mkdir -p $RPM_BUILD_ROOT/opt/DSR
dos2unix *
cp -R * $RPM_BUILD_ROOT/opt/DSR
mkdir -p $RPM_BUILD_ROOT/etc/profile.d
cp setup/dsr.sh $RPM_BUILD_ROOT/etc/profile.d


%files
%doc /opt/DSR/manuals/DSR-manual.pdf
#%defattr(644, root, root, 755)
%attr(644, root, root) /etc/profile.d/dsr.sh
%config /etc/profile.d/dsr.sh
#%attr(666, root, users) /opt/DSR/dsr_user_db.txt
%attr(755, root, users) /opt/DSR/dsr
/opt/DSR/afix.py
/opt/DSR/atomhandling.py
/opt/DSR/atoms.py
/opt/DSR/setup/dsr.sh
/opt/DSR/constants.py
/opt/DSR/dbfile.py
/opt/DSR/dsr.py
/opt/DSR/dsrparse.py
/opt/DSR/export.py
/opt/DSR/misc.py
/opt/DSR/terminalsize.py
/opt/DSR/options.py
/opt/DSR/resfile.py
/opt/DSR/restraints.py
/opt/DSR/resi.py
/opt/DSR/refine.py
/opt/DSR/pyperclip.py
/opt/DSR/dsr_db.txt
#/opt/DSR/dsr_user_db.txt
/opt/DSR/manuals/DSR-manual.pdf
#/opt/DSR/setup/set_environ.py
/opt/DSR/setup/OlexDSR.py
/opt/DSR/setup/custom.xld
/opt/DSR/networkx/
%attr(777, root, users) /opt/DSR/example
%attr(666, root, users) /opt/DSR/example/p21c.hkl
%attr(666, root, users) /opt/DSR/example/p21c.res
%attr(666, root, users) /opt/DSR/example/p21c_step0.res
%attr(666, root, users) /opt/DSR/example/p21c_step1.res
%attr(666, root, users) /opt/DSR/example/p21c_step2.res
%attr(666, root, users) /opt/DSR/example/p21c_step3.res
%attr(666, root, users) /opt/DSR/example/p21c-step2.ins
