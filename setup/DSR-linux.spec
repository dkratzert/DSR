Summary: DSR - A program for modelling of disordered solvents with SHELXL
Name: DSR
Version: 217
Release: 0
BuildArch: noarch
URL: https://www.xs3.uni-freiburg.de/research/dsr
License: Beerware
Group: Sciences/Chemistry
Source: %name-%version.tar.gz

PreReq: xclip %py_requires python2-numpy

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


%build
cat > %name.sh << EOF
#!/bin/sh
export DSR_DIR=%_datadir/DSR
PYTHON_EXE=\$(which python)
if [ \$# -eq 0 ]; then
    \$PYTHON_EXE \$DSR_DIR/dsr.py --help
else
    \$PYTHON_EXE \$DSR_DIR/dsr.py \$*
fi
EOF

%install
mkdir -p %buildroot%_bindir
mkdir -p %buildroot%_datadir/%name
mkdir -p %buildroot%_datadir/%name/manuals
mkdir -p %buildroot%_datadir/%name/example
mkdir -p %buildroot%_datadir/%name/mpmath
mkdir -p %buildroot%_datadir/%name/mpmath/calculus
mkdir -p %buildroot%_datadir/%name/mpmath/functions
mkdir -p %buildroot%_datadir/%name/mpmath/libmp
mkdir -p %buildroot%_datadir/%name/mpmath/matrices
mkdir -p %buildroot%_datadir/%name/rmsd
mkdir -p %buildroot%_datadir/%name/networkx
mkdir -p %buildroot%_datadir/%name/networkx/external/decorator/decorator3
mkdir -p %buildroot%_datadir/%name/networkx/external/decorator/decorator2
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/assortativity/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/link_analysis/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/shortest_paths/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/approximation/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/chordal/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/components/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/centrality/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/flow/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/community/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/connectivity/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/bipartite/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/traversal/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/operators/tests
mkdir -p %buildroot%_datadir/%name/networkx/algorithms/isomorphism/tests
mkdir -p %buildroot%_datadir/%name/networkx/tests
mkdir -p %buildroot%_datadir/%name/networkx/generators/tests
mkdir -p %buildroot%_datadir/%name/networkx/linalg/tests
mkdir -p %buildroot%_datadir/%name/networkx/readwrite/json_graph/tests
mkdir -p %buildroot%_datadir/%name/networkx/readwrite/tests
mkdir -p %buildroot%_datadir/%name/networkx/classes/tests
mkdir -p %buildroot%_datadir/%name/networkx/utils/tests
mkdir -p %buildroot%_datadir/%name/networkx/drawing/tests
cp -R networkx %buildroot%_datadir/%name/networkx
cp -R rmsd %buildroot%_datadir/%name/rmsd
cp -R mpmath %buildroot%_datadir/%name/mpmath

install -m 755 %name.sh %buildroot%_bindir/dsr
install -m 644 *.py %buildroot%_datadir/%name
install -m 644 dsr_db.txt %buildroot%_datadir/%name
install -m 644 manuals/DSR-manual.pdf %buildroot%_datadir/%name/manuals
install -m 644 example/* %buildroot%_datadir/%name/example
install -mr 644 mpmath/* %buildroot%_datadir/%name/mpmath
install -m 644 rmsd/* %buildroot%_datadir/%name/rmsd
install -m 644 networkx/* %buildroot%_datadir/%name/networkx

%files
%doc README
%_bindir/dsr
%_datadir/%name
%_datadir/rmsd/*
%_datadir/mpmath/*
%_datadir/networkx/*

%changelog
* Sun Aug 19 2018 Daniel Kratzert
* Thu Jan 04 2018 Denis G. Samsonenko <ogion@altlinux.org> 205-alt1
- initial build for ALT
