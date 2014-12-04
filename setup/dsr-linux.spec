Summary: DSR - A program for modelling of disordered solvents with SHELXL
Name: DSR
Provides: DSR
Packager: dkratzert@gmx.de
Version: 1.5.13
Release: 0
Requires: python, xclip
Prefix: /opt
BuildRoot: %{_tmppath}/%{name}-%{version}
BuildArch: noarch
URL: https://www.xs3.uni-freiburg.de/research/dsr
License: Beerware
Source: DSR-%{version}.tar.gz
Group: Productivity/Scientific/Chemistry

%define _unpackaged_files_terminate_build 0

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
rm -fr $RPM_BUILD_ROOT
mkdir %{buildroot}
mkdir %{buildroot}/DSR
mkdir %{buildroot}/DSR/manuals
mkdir %{buildroot}/DSR/example
mkdir %{buildroot}/DSR/setup
mkdir -p %{buildroot}/etc/profile.d
touch %{buildroot}/DSR/dsr_user_db.txt

%post
chmod a+rw %{prefix}/DSR/example
source /etc/profile

%install
install -m 644 afix.py %{buildroot}/DSR/afix.py
install -m 644 atomhandling.py %{buildroot}/DSR/atomhandling.py
install -m 644 atoms.py %{buildroot}/DSR/atoms.py
install -m 644 constants.py %{buildroot}/DSR/constants.py
install -m 644 dbfile.py %{buildroot}/DSR/dbfile.py
install -m 644 dsr.py %{buildroot}/DSR/dsr.py
install -m 755 dsr %{buildroot}/DSR/dsr
install -m 644 dsrparse.py %{buildroot}/DSR/dsrparse.py
install -m 644 export.py %{buildroot}/DSR/export.py
install -m 644 misc.py %{buildroot}/DSR/misc.py
install -m 644 options.py %{buildroot}/DSR/options.py
install -m 644 terminalsize.py %{buildroot}/DSR/terminalsize.py
install -m 644 resfile.py %{buildroot}/DSR/resfile.py
install -m 644 restraints.py %{buildroot}/DSR/restraints.py
install -m 644 resi.py %{buildroot}/DSR/resi.py
install -m 644 refine.py %{buildroot}/DSR/refine.py
install -m 644 pyperclip.py %{buildroot}/DSR/pyperclip.py
install -m 644 dsr_db.txt %{buildroot}/DSR/dsr_db.txt
install -m 644 manuals/DSR-manual.pdf %{buildroot}/DSR/manuals/DSR-manual.pdf
install -m 644 setup/dsr.sh %{buildroot}/etc/profile.d/dsr.sh
install -m 666 example/p21c.hkl %{buildroot}/DSR/example/p21c.hkl
install -m 666 example/p21c.res %{buildroot}/DSR/example/p21c.res
install -m 666 example/p21c_step0.res %{buildroot}/DSR/example/p21c_step0.res
install -m 666 example/p21c_step1.res %{buildroot}/DSR/example/p21c_step1.res
install -m 666 example/p21c_step2.res %{buildroot}/DSR/example/p21c_step2.res
install -m 666 example/p21c_step3.res %{buildroot}/DSR/example/p21c_step3.res
install -m 666 example/p21c-step2.ins %{buildroot}/DSR/example/p21c-step2.ins
install -d -m 755 networkx/ %{buildroot}/DSR/networkx
dos2unix %{buildroot}/*

%files
#/DSR/networkx/*
%doc /DSR/manuals/DSR-manual.pdf
/DSR/dsr
/DSR/afix.py
/DSR/atomhandling.py
/DSR/atoms.py
#/DSR/setup/dsr.sh
%config /etc/profile.d/dsr.sh
/DSR/constants.py
/DSR/dbfile.py
/DSR/dsr.py
/DSR/dsr
/DSR/dsrparse.py
/DSR/export.py
/DSR/misc.py
/DSR/terminalsize.py
/DSR/options.py
/DSR/resfile.py 
/DSR/restraints.py
/DSR/resi.py
/DSR/refine.py
/DSR/pyperclip.py
/DSR/dsr_db.txt
/DSR/dsr_user_db.txt
/DSR/manuals/DSR-manual.pdf
#/DSR/setup/set_environ.py
/DSR/networkx/
/DSR/example/
/DSR/example/p21c.hkl
/DSR/example/p21c.res
/DSR/example/p21c_step0.res
/DSR/example/p21c_step1.res
/DSR/example/p21c_step2.res
/DSR/example/p21c_step3.res
/DSR/example/p21c-step2.ins
