Summary: DSR - A program for modelling of disordered solvents with SHELXL
Name: DSR
Provides: DSR
Packager: dkratzert@gmx.de
Version: 1.5.13
Release: 0
Requires: python, xclip
Prefix: /opt
BuildRoot: %{_tmppath}/%{name}-%{version}-build
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
mkdir -p %{prefix}/DSR
mkdir -p %{prefix}/DSR/manuals
mkdir -p %{prefix}/DSR/example
mkdir -p %{prefix}/DSR/setup

%post
if [ ! -f %{prefix}/DSR/dsr_user_db.txt ]
then
    touch %{prefix}/DSR/dsr_user_db.txt
    chmod a+rw %{prefix}/DSR/dsr_user_db.txt
else
    chmod a+rw %{prefix}/DSR/dsr_user_db.txt
fi
chmod a+rw %{prefix}/DSR/example
source /etc/profile

%install
install -m 644 afix.py %{prefix}/DSR/afix.py
install -m 644 atomhandling.py %{prefix}/DSR/atomhandling.py
install -m 644 atoms.py %{prefix}/DSR/atoms.py
install -m 644 constants.py %{prefix}/DSR/constants.py
install -m 644 dbfile.py %{prefix}/DSR/dbfile.py
install -m 644 dsr.py %{prefix}/DSR/dsr.py
install -m 755 dsr %{prefix}/DSR/dsr
install -m 644 dsrparse.py %{prefix}/DSR/dsrparse.py
install -m 644 export.py %{prefix}/DSR/export.py
install -m 644 misc.py %{prefix}/DSR/misc.py
install -m 644 options.py %{prefix}/DSR/options.py
install -m 644 terminalsize.py %{prefix}/DSR/terminalsize.py
install -m 644 resfile.py %{prefix}/DSR/resfile.py
install -m 644 restraints.py %{prefix}/DSR/restraints.py
install -m 644 resi.py %{prefix}/DSR/resi.py
install -m 644 refine.py %{prefix}/DSR/refine.py
install -m 644 pyperclip.py %{prefix}/DSR/pyperclip.py
install -m 644 dsr_db.txt %{prefix}/DSR/dsr_db.txt
install -m 644 manuals/DSR-manual.pdf %{prefix}/DSR/manuals/DSR-manual.pdf
install -m 644 setup/dsr.sh /etc/profile.d/dsr.sh
install -m 644 example/p21c.hkl %{prefix}/DSR/example/p21c.hkl
install -m 644 example/p21c.res %{prefix}/DSR/example/p21c.res
install -m 644 example/p21c_step0.res %{prefix}/DSR/example/p21c_step0.res
install -m 644 example/p21c_step1.res %{prefix}/DSR/example/p21c_step1.res
install -m 644 example/p21c_step2.res %{prefix}/DSR/example/p21c_step2.res
install -m 644 example/p21c_step3.res %{prefix}/DSR/example/p21c_step3.res
install -m 644 example/p21c-step2.ins %{prefix}/DSR/example/p21c-step2.ins
install -d -m 755 networkx %{prefix}/DSR/networkx
mkdir -p $RPM_BUILD_ROOT%{prefix}/DSR
dos2unix *
cp -R * $RPM_BUILD_ROOT%{prefix}/DSR
mkdir -p $RPM_BUILD_ROOT/etc/profile.d
cp setup/dsr.sh $RPM_BUILD_ROOT/etc/profile.d


%files
%doc %{prefix}/DSR/manuals/DSR-manual.pdf
#%defattr(644, root, root, 755)
%attr(644, root, root) /etc/profile.d/dsr.sh
%config /etc/profile.d/dsr.sh
#%attr(666, root, users) %{prefix}/DSR/dsr_user_db.txt
%attr(755, root, users) %{prefix}/DSR/dsr
%{prefix}/DSR/afix.py
%{prefix}/DSR/atomhandling.py
%{prefix}/DSR/atoms.py
%{prefix}/DSR/setup/dsr.sh
%{prefix}/DSR/constants.py
%{prefix}/DSR/dbfile.py
%{prefix}/DSR/dsr.py
%{prefix}/DSR/dsrparse.py
%{prefix}/DSR/export.py
%{prefix}/DSR/misc.py
%{prefix}/DSR/terminalsize.py
%{prefix}/DSR/options.py
%{prefix}/DSR/resfile.py
%{prefix}/DSR/restraints.py
%{prefix}/DSR/resi.py
%{prefix}/DSR/refine.py
%{prefix}/DSR/pyperclip.py
%{prefix}/DSR/dsr_db.txt
#%{prefix}/DSR/dsr_user_db.txt
%{prefix}/DSR/manuals/DSR-manual.pdf
#%{prefix}/DSR/setup/set_environ.py
%{prefix}/DSR/networkx/
%attr(777, root, users) %{prefix}/DSR/example
%attr(666, root, users) %{prefix}/DSR/example/p21c.hkl
%attr(666, root, users) %{prefix}/DSR/example/p21c.res
%attr(666, root, users) %{prefix}/DSR/example/p21c_step0.res
%attr(666, root, users) %{prefix}/DSR/example/p21c_step1.res
%attr(666, root, users) %{prefix}/DSR/example/p21c_step2.res
%attr(666, root, users) %{prefix}/DSR/example/p21c_step3.res
%attr(666, root, users) %{prefix}/DSR/example/p21c-step2.ins
