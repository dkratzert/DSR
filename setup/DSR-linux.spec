Summary: DSR - A program for modelling of disordered solvents with SHELXL
Name: DSR
Provides: DSR
Packager: dkratzert@gmx.de
Version: 234
Release: 0
Requires: python
Prefix: /opt
BuildRoot: %{_tmppath}/%{name}-%{version}
BuildArch: noarch
URL: https://dkratzert.de/dsr.html
License: Beerware
Source: %{name}-%{version}.tar.gz
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
mkdir -p %{buildroot}%{prefix}/DSR/networkx
mkdir -p %{buildroot}%{prefix}/DSR/manuals
mkdir -p %{buildroot}%{prefix}/DSR/example
mkdir -p %{buildroot}%{prefix}/DSR/setup
mkdir -p %{buildroot}%{prefix}/DSR/fit
mkdir -p %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator3
mkdir -p %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator2
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/link_analysis
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/chordal
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/components
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/flow
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/community
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/connectivity
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/traversal
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/operators
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism
mkdir -p %{buildroot}%{prefix}/DSR/networkx
mkdir -p %{buildroot}%{prefix}/DSR/networkx/generators
mkdir -p %{buildroot}%{prefix}/DSR/networkx/linalg
mkdir -p %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph
mkdir -p %{buildroot}%{prefix}/DSR/networkx/readwrite
mkdir -p %{buildroot}%{prefix}/DSR/networkx/classes
mkdir -p %{buildroot}%{prefix}/DSR/networkx/utils
mkdir -p %{buildroot}%{prefix}/DSR/networkx/drawing
mkdir -p %{buildroot}%{prefix}/DSR/networkx/testing
mkdir -p %{buildroot}%{prefix}/DSR/mpmath/calculus
mkdir -p %{buildroot}%{prefix}/DSR/mpmath/functions
mkdir -p %{buildroot}%{prefix}/DSR/mpmath/libmp
mkdir -p %{buildroot}%{prefix}/DSR/mpmath/matrices

mkdir -p %{buildroot}/etc/profile.d
chmod a+rw %{buildroot}%{prefix}/DSR/example

%post
echo "DSR was installed in /opt/DSR"
echo "Please restart your system to ensure that the environment variables are set."
source /etc/profile.d/dsr.sh



%install
install -m 644 README %{buildroot}%{prefix}/DSR/README
install -m 644 afix.py %{buildroot}%{prefix}/DSR/afix.py
install -m 644 atomhandling.py %{buildroot}%{prefix}/DSR/atomhandling.py
install -m 644 atoms.py %{buildroot}%{prefix}/DSR/atoms.py
install -m 644 constants.py %{buildroot}%{prefix}/DSR/constants.py
install -m 644 dbfile.py %{buildroot}%{prefix}/DSR/dbfile.py
install -m 644 dsr.py %{buildroot}%{prefix}/DSR/dsr.py
install -m 755 dsr %{buildroot}%{prefix}/DSR/dsr
install -m 644 dsrparse.py %{buildroot}%{prefix}/DSR/dsrparse.py
install -m 644 export.py %{buildroot}%{prefix}/DSR/export.py
install -m 644 misc.py %{buildroot}%{prefix}/DSR/misc.py
install -m 644 options.py %{buildroot}%{prefix}/DSR/options.py
install -m 644 terminalsize.py %{buildroot}%{prefix}/DSR/terminalsize.py
install -m 644 resfile.py %{buildroot}%{prefix}/DSR/resfile.py
install -m 644 cf3fit.py %{buildroot}%{prefix}/DSR/cf3fit.py
install -m 644 selfupdate.py %{buildroot}%{prefix}/DSR/selfupdate.py
install -m 644 elements.py %{buildroot}%{prefix}/DSR/elements.py
install -m 644 restraints.py %{buildroot}%{prefix}/DSR/restraints.py
install -m 644 resi.py %{buildroot}%{prefix}/DSR/resi.py
install -m 644 refine.py %{buildroot}%{prefix}/DSR/refine.py
install -m 644 pyperclip.py %{buildroot}%{prefix}/DSR/pyperclip.py
install -m 644 dsr_db.txt %{buildroot}%{prefix}/DSR/dsr_db.txt
install -m 644 manuals/DSR-manual.pdf %{buildroot}%{prefix}/DSR/manuals/DSR-manual.pdf
install -m 644 setup/dsr.sh %{buildroot}/etc/profile.d/dsr.sh
install -m 666 example/p21c.hkl %{buildroot}%{prefix}/DSR/example/p21c.hkl
install -m 666 example/p21c.res %{buildroot}%{prefix}/DSR/example/p21c.res
install -m 666 example/p21c_step0.res %{buildroot}%{prefix}/DSR/example/p21c_step0.res
install -m 666 example/p21c_step1.res %{buildroot}%{prefix}/DSR/example/p21c_step1.res
install -m 666 example/p21c_step2.ins %{buildroot}%{prefix}/DSR/example/p21c_step2.ins
install -m 666 example/p21c_step3.res %{buildroot}%{prefix}/DSR/example/p21c_step3.res
install -m 666 example/p21c_final.res %{buildroot}%{prefix}/DSR/example/p21c_final.res
install -m 666 example/p21n_cf3.hkl %{buildroot}%{prefix}/DSR/example/p21n_cf3.hkl
install -m 666 example/p21n_cf3.res %{buildroot}%{prefix}/DSR/example/p21n_cf3.res

install -m 644 fit/__init__.py %{buildroot}%{prefix}/DSR/fit/__init__.py
install -m 644 fit/quatfit.py %{buildroot}%{prefix}/DSR/fit/quatfit.py
install -m 644 fit/LICENSE.txt %{buildroot}%{prefix}/DSR/fit/LICENSE.txt

install -m 644 networkx/external/decorator/decorator3/_decorator3.py            %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator3/_decorator3.py
install -m 644 networkx/external/decorator/decorator3/__init__.py               %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator3/__init__.py
install -m 644 networkx/external/decorator/__init__.py                          %{buildroot}%{prefix}/DSR/networkx/external/decorator/__init__.py
install -m 644 networkx/external/decorator/decorator2/_decorator2.py            %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator2/_decorator2.py
install -m 644 networkx/external/decorator/decorator2/__init__.py               %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator2/__init__.py
install -m 644 networkx/external/__init__.py                                    %{buildroot}%{prefix}/DSR/networkx/external/__init__.py
install -m 644 networkx/algorithms/richclub.py                                  %{buildroot}%{prefix}/DSR/networkx/algorithms/richclub.py
install -m 644 networkx/algorithms/assortativity/pairs.py                       %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/pairs.py
install -m 644 networkx/algorithms/assortativity/correlation.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/correlation.py
install -m 644 networkx/algorithms/assortativity/mixing.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/mixing.py
install -m 644 networkx/algorithms/assortativity/__init__.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/__init__.py
install -m 644 networkx/algorithms/assortativity/neighbor_degree.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/neighbor_degree.py
install -m 644 networkx/algorithms/assortativity/connectivity.py                %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/connectivity.py
install -m 644 networkx/algorithms/link_analysis/pagerank_alg.py                %{buildroot}%{prefix}/DSR/networkx/algorithms/link_analysis/pagerank_alg.py
install -m 644 networkx/algorithms/link_analysis/__init__.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/link_analysis/__init__.py
install -m 644 networkx/algorithms/link_analysis/hits_alg.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/link_analysis/hits_alg.py
install -m 644 networkx/algorithms/smetric.py                                   %{buildroot}%{prefix}/DSR/networkx/algorithms/smetric.py
install -m 644 networkx/algorithms/graphical.py                                 %{buildroot}%{prefix}/DSR/networkx/algorithms/graphical.py
install -m 644 networkx/algorithms/shortest_paths/dense.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/dense.py
install -m 644 networkx/algorithms/shortest_paths/unweighted.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/unweighted.py
install -m 644 networkx/algorithms/shortest_paths/weighted.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/weighted.py
install -m 644 networkx/algorithms/shortest_paths/astar.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/astar.py
install -m 644 networkx/algorithms/shortest_paths/__init__.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/__init__.py
install -m 644 networkx/algorithms/shortest_paths/generic.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/generic.py
install -m 644 networkx/algorithms/simple_paths.py                              %{buildroot}%{prefix}/DSR/networkx/algorithms/simple_paths.py
install -m 644 networkx/algorithms/vitality.py                                  %{buildroot}%{prefix}/DSR/networkx/algorithms/vitality.py
install -m 644 networkx/algorithms/approximation/independent_set.py               %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/independent_set.py
install -m 644 networkx/algorithms/approximation/dominating_set.py                %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/dominating_set.py
install -m 644 networkx/algorithms/approximation/matching.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/matching.py
install -m 644 networkx/algorithms/approximation/ramsey.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/ramsey.py
install -m 644 networkx/algorithms/approximation/vertex_cover.py                  %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/vertex_cover.py
install -m 644 networkx/algorithms/approximation/clique.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/clique.py
install -m 644 networkx/algorithms/approximation/__init__.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/__init__.py
install -m 644 networkx/algorithms/chordal/chordal_alg.py                         %{buildroot}%{prefix}/DSR/networkx/algorithms/chordal/chordal_alg.py
install -m 644 networkx/algorithms/chordal/__init__.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/chordal/__init__.py
install -m 644 networkx/algorithms/hierarchy.py                                   %{buildroot}%{prefix}/DSR/networkx/algorithms/hierarchy.py
install -m 644 networkx/algorithms/block.py                                       %{buildroot}%{prefix}/DSR/networkx/algorithms/block.py
install -m 644 networkx/algorithms/core.py                                        %{buildroot}%{prefix}/DSR/networkx/algorithms/core.py
install -m 644 networkx/algorithms/distance_regular.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/distance_regular.py
install -m 644 networkx/algorithms/components/strongly_connected.py               %{buildroot}%{prefix}/DSR/networkx/algorithms/components/strongly_connected.py
install -m 644 networkx/algorithms/components/weakly_connected.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/components/weakly_connected.py
install -m 644 networkx/algorithms/components/attracting.py                       %{buildroot}%{prefix}/DSR/networkx/algorithms/components/attracting.py
install -m 644 networkx/algorithms/components/biconnected.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/components/biconnected.py
install -m 644 networkx/algorithms/components/__init__.py                         %{buildroot}%{prefix}/DSR/networkx/algorithms/components/__init__.py
install -m 644 networkx/algorithms/components/connected.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/components/connected.py
install -m 644 networkx/algorithms/centrality/closeness.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/closeness.py
install -m 644 networkx/algorithms/centrality/katz.py                                  %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/katz.py
install -m 644 networkx/algorithms/centrality/betweenness.py                           %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/betweenness.py
install -m 644 networkx/algorithms/centrality/flow_matrix.py                           %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/flow_matrix.py
install -m 644 networkx/algorithms/centrality/load.py                                  %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/load.py
install -m 644 networkx/algorithms/centrality/current_flow_closeness.py                %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/current_flow_closeness.py
install -m 644 networkx/algorithms/centrality/communicability_alg.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/communicability_alg.py
install -m 644 networkx/algorithms/centrality/eigenvector.py                           %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/eigenvector.py
install -m 644 networkx/algorithms/centrality/degree_alg.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/degree_alg.py
install -m 644 networkx/algorithms/centrality/current_flow_betweenness_subset.py       %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/current_flow_betweenness_subset.py
install -m 644 networkx/algorithms/centrality/current_flow_betweenness.py              %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/current_flow_betweenness.py
install -m 644 networkx/algorithms/centrality/betweenness_subset.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/betweenness_subset.py
install -m 644 networkx/algorithms/centrality/__init__.py                              %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/__init__.py
install -m 644 networkx/algorithms/matching.py                                         %{buildroot}%{prefix}/DSR/networkx/algorithms/matching.py
install -m 644 networkx/algorithms/dag.py                                              %{buildroot}%{prefix}/DSR/networkx/algorithms/dag.py
install -m 644 networkx/algorithms/boundary.py                                         %{buildroot}%{prefix}/DSR/networkx/algorithms/boundary.py
install -m 644 networkx/algorithms/euler.py                                            %{buildroot}%{prefix}/DSR/networkx/algorithms/euler.py
install -m 644 networkx/algorithms/flow/mincost.py                                     %{buildroot}%{prefix}/DSR/networkx/algorithms/flow/mincost.py
install -m 644 networkx/algorithms/flow/maxflow.py                                     %{buildroot}%{prefix}/DSR/networkx/algorithms/flow/maxflow.py
install -m 644 networkx/algorithms/flow/__init__.py                                    %{buildroot}%{prefix}/DSR/networkx/algorithms/flow/__init__.py
install -m 644 networkx/algorithms/distance_measures.py                                %{buildroot}%{prefix}/DSR/networkx/algorithms/distance_measures.py
install -m 644 networkx/algorithms/isolate.py                                          %{buildroot}%{prefix}/DSR/networkx/algorithms/isolate.py
install -m 644 networkx/algorithms/community/kclique.py                                %{buildroot}%{prefix}/DSR/networkx/algorithms/community/kclique.py
install -m 644 networkx/algorithms/community/__init__.py                               %{buildroot}%{prefix}/DSR/networkx/algorithms/community/__init__.py
install -m 644 networkx/algorithms/swap.py                                             %{buildroot}%{prefix}/DSR/networkx/algorithms/swap.py
install -m 644 networkx/algorithms/cycles.py                                           %{buildroot}%{prefix}/DSR/networkx/algorithms/cycles.py
install -m 644 networkx/algorithms/connectivity/__init__.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/connectivity/__init__.py
install -m 644 networkx/algorithms/connectivity/connectivity.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/connectivity/connectivity.py
install -m 644 networkx/algorithms/connectivity/cuts.py                                %{buildroot}%{prefix}/DSR/networkx/algorithms/connectivity/cuts.py
install -m 644 networkx/algorithms/clique.py                                           %{buildroot}%{prefix}/DSR/networkx/algorithms/clique.py
install -m 644 networkx/algorithms/__init__.py                                         %{buildroot}%{prefix}/DSR/networkx/algorithms/__init__.py
install -m 644 networkx/algorithms/mst.py                                              %{buildroot}%{prefix}/DSR/networkx/algorithms/mst.py
install -m 644 networkx/algorithms/bipartite/centrality.py                             %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/centrality.py
install -m 644 networkx/algorithms/bipartite/projection.py                             %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/projection.py
install -m 644 networkx/algorithms/bipartite/basic.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/basic.py
install -m 644 networkx/algorithms/bipartite/spectral.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/spectral.py
install -m 644 networkx/algorithms/bipartite/__init__.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/__init__.py
install -m 644 networkx/algorithms/bipartite/redundancy.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/redundancy.py
install -m 644 networkx/algorithms/bipartite/cluster.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/cluster.py
install -m 644 networkx/algorithms/traversal/depth_first_search.py         %{buildroot}%{prefix}/DSR/networkx/algorithms/traversal/depth_first_search.py
install -m 644 networkx/algorithms/traversal/__init__.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/traversal/__init__.py
install -m 644 networkx/algorithms/traversal/breadth_first_search.py       %{buildroot}%{prefix}/DSR/networkx/algorithms/traversal/breadth_first_search.py
install -m 644 networkx/algorithms/mis.py                                  %{buildroot}%{prefix}/DSR/networkx/algorithms/mis.py
install -m 644 networkx/algorithms/operators/unary.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/unary.py
install -m 644 networkx/algorithms/operators/all.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/all.py
install -m 644 networkx/algorithms/operators/product.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/product.py
install -m 644 networkx/algorithms/operators/binary.py                     %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/binary.py
install -m 644 networkx/algorithms/operators/__init__.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/__init__.py
install -m 644 networkx/algorithms/isomorphism/isomorphvf2.py              %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/isomorphvf2.py
install -m 644 networkx/algorithms/isomorphism/matchhelpers.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/matchhelpers.py
install -m 644 networkx/algorithms/isomorphism/__init__.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/__init__.py
install -m 644 networkx/algorithms/isomorphism/isomorph.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/isomorph.py
install -m 644 networkx/algorithms/isomorphism/vf2userfunc.py              %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/vf2userfunc.py
install -m 644 networkx/algorithms/cluster.py                              %{buildroot}%{prefix}/DSR/networkx/algorithms/cluster.py
install -m 644 networkx/generators/classic.py                              %{buildroot}%{prefix}/DSR/networkx/generators/classic.py
install -m 644 networkx/generators/social.py                               %{buildroot}%{prefix}/DSR/networkx/generators/social.py
install -m 644 networkx/generators/intersection.py                         %{buildroot}%{prefix}/DSR/networkx/generators/intersection.py
install -m 644 networkx/generators/threshold.py                            %{buildroot}%{prefix}/DSR/networkx/generators/threshold.py
install -m 644 networkx/generators/ego.py                                  %{buildroot}%{prefix}/DSR/networkx/generators/ego.py
install -m 644 networkx/generators/atlas.py                                %{buildroot}%{prefix}/DSR/networkx/generators/atlas.py
install -m 644 networkx/generators/hybrid.py                               %{buildroot}%{prefix}/DSR/networkx/generators/hybrid.py
install -m 644 networkx/generators/directed.py                             %{buildroot}%{prefix}/DSR/networkx/generators/directed.py
install -m 644 networkx/generators/random_graphs.py                        %{buildroot}%{prefix}/DSR/networkx/generators/random_graphs.py
install -m 644 networkx/generators/bipartite.py                            %{buildroot}%{prefix}/DSR/networkx/generators/bipartite.py
install -m 644 networkx/generators/line.py                                 %{buildroot}%{prefix}/DSR/networkx/generators/line.py
install -m 644 networkx/generators/__init__.py                             %{buildroot}%{prefix}/DSR/networkx/generators/__init__.py
install -m 644 networkx/generators/geometric.py                            %{buildroot}%{prefix}/DSR/networkx/generators/geometric.py
install -m 644 networkx/generators/random_clustered.py                     %{buildroot}%{prefix}/DSR/networkx/generators/random_clustered.py
install -m 644 networkx/generators/small.py                                %{buildroot}%{prefix}/DSR/networkx/generators/small.py
install -m 644 networkx/generators/degree_seq.py                           %{buildroot}%{prefix}/DSR/networkx/generators/degree_seq.py
install -m 644 networkx/generators/stochastic.py                           %{buildroot}%{prefix}/DSR/networkx/generators/stochastic.py
install -m 644 networkx/exception.py                                       %{buildroot}%{prefix}/DSR/networkx/exception.py
install -m 644 networkx/version.py                                         %{buildroot}%{prefix}/DSR/networkx/version.py
install -m 644 networkx/release.py                                         %{buildroot}%{prefix}/DSR/networkx/release.py
install -m 644 networkx/linalg/attrmatrix.py                               %{buildroot}%{prefix}/DSR/networkx/linalg/attrmatrix.py
install -m 644 networkx/linalg/graphmatrix.py                              %{buildroot}%{prefix}/DSR/networkx/linalg/graphmatrix.py
install -m 644 networkx/linalg/__init__.py                                 %{buildroot}%{prefix}/DSR/networkx/linalg/__init__.py
install -m 644 networkx/linalg/laplacianmatrix.py                          %{buildroot}%{prefix}/DSR/networkx/linalg/laplacianmatrix.py
install -m 644 networkx/linalg/spectrum.py                                 %{buildroot}%{prefix}/DSR/networkx/linalg/spectrum.py
install -m 644 networkx/readwrite/leda.py                                  %{buildroot}%{prefix}/DSR/networkx/readwrite/leda.py
install -m 644 networkx/readwrite/multiline_adjlist.py                     %{buildroot}%{prefix}/DSR/networkx/readwrite/multiline_adjlist.py
install -m 644 networkx/readwrite/json_graph/tree.py                       %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph/tree.py
install -m 644 networkx/readwrite/json_graph/serialize.py                  %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph/serialize.py
install -m 644 networkx/readwrite/json_graph/node_link.py                  %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph/node_link.py
install -m 644 networkx/readwrite/json_graph/__init__.py                   %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph/__init__.py
install -m 644 networkx/readwrite/json_graph/adjacency.py                  %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph/adjacency.py
install -m 644 networkx/readwrite/edgelist.py                              %{buildroot}%{prefix}/DSR/networkx/readwrite/edgelist.py
install -m 644 networkx/readwrite/gml.py                                   %{buildroot}%{prefix}/DSR/networkx/readwrite/gml.py
install -m 644 networkx/readwrite/pajek.py                                 %{buildroot}%{prefix}/DSR/networkx/readwrite/pajek.py
install -m 644 networkx/readwrite/graphml.py                               %{buildroot}%{prefix}/DSR/networkx/readwrite/graphml.py
install -m 644 networkx/readwrite/gexf.py                                  %{buildroot}%{prefix}/DSR/networkx/readwrite/gexf.py
install -m 644 networkx/readwrite/adjlist.py                               %{buildroot}%{prefix}/DSR/networkx/readwrite/adjlist.py
install -m 644 networkx/readwrite/gpickle.py                               %{buildroot}%{prefix}/DSR/networkx/readwrite/gpickle.py
install -m 644 networkx/readwrite/__init__.py                              %{buildroot}%{prefix}/DSR/networkx/readwrite/__init__.py
install -m 644 networkx/readwrite/sparsegraph6.py                          %{buildroot}%{prefix}/DSR/networkx/readwrite/sparsegraph6.py
install -m 644 networkx/readwrite/nx_yaml.py                               %{buildroot}%{prefix}/DSR/networkx/readwrite/nx_yaml.py
install -m 644 networkx/readwrite/p2g.py                                   %{buildroot}%{prefix}/DSR/networkx/readwrite/p2g.py
install -m 644 networkx/readwrite/nx_shp.py                                %{buildroot}%{prefix}/DSR/networkx/readwrite/nx_shp.py
install -m 644 networkx/relabel.py                                         %{buildroot}%{prefix}/DSR/networkx/relabel.py
install -m 644 networkx/classes/graph.py                                   %{buildroot}%{prefix}/DSR/networkx/classes/graph.py
install -m 644 networkx/classes/multigraph.py                              %{buildroot}%{prefix}/DSR/networkx/classes/multigraph.py
install -m 644 networkx/classes/function.py                                %{buildroot}%{prefix}/DSR/networkx/classes/function.py
install -m 644 networkx/classes/digraph.py                                 %{buildroot}%{prefix}/DSR/networkx/classes/digraph.py
install -m 644 networkx/classes/__init__.py                                %{buildroot}%{prefix}/DSR/networkx/classes/__init__.py
install -m 644 networkx/classes/multidigraph.py                            %{buildroot}%{prefix}/DSR/networkx/classes/multidigraph.py
install -m 644 networkx/utils/union_find.py                                %{buildroot}%{prefix}/DSR/networkx/utils/union_find.py
install -m 644 networkx/utils/decorators.py                                %{buildroot}%{prefix}/DSR/networkx/utils/decorators.py
install -m 644 networkx/utils/rcm.py                                       %{buildroot}%{prefix}/DSR/networkx/utils/rcm.py
install -m 644 networkx/utils/misc.py                                      %{buildroot}%{prefix}/DSR/networkx/utils/misc.py
install -m 644 networkx/utils/random_sequence.py                           %{buildroot}%{prefix}/DSR/networkx/utils/random_sequence.py
install -m 644 networkx/utils/__init__.py                                  %{buildroot}%{prefix}/DSR/networkx/utils/__init__.py
install -m 644 networkx/drawing/nx_agraph.py                               %{buildroot}%{prefix}/DSR/networkx/drawing/nx_agraph.py
install -m 644 networkx/drawing/nx_pydot.py                                %{buildroot}%{prefix}/DSR/networkx/drawing/nx_pydot.py
install -m 644 networkx/drawing/__init__.py                                %{buildroot}%{prefix}/DSR/networkx/drawing/__init__.py
install -m 644 networkx/drawing/layout.py                                  %{buildroot}%{prefix}/DSR/networkx/drawing/layout.py
install -m 644 networkx/drawing/nx_pylab.py                                %{buildroot}%{prefix}/DSR/networkx/drawing/nx_pylab.py
install -m 644 networkx/__init__.py                                        %{buildroot}%{prefix}/DSR/networkx/__init__.py
install -m 644 networkx/convert.py                                         %{buildroot}%{prefix}/DSR/networkx/convert.py
install -m 644 networkx/testing/utils.py                                   %{buildroot}%{prefix}/DSR/networkx/testing/utils.py
install -m 644 networkx/testing/__init__.py                                %{buildroot}%{prefix}/DSR/networkx/testing/__init__.py

install -m 644 mpmath/conftest.py                   %{buildroot}%{prefix}/DSR/mpmath/conftest.py
install -m 644 mpmath/ctx_base.py                   %{buildroot}%{prefix}/DSR/mpmath/ctx_base.py
install -m 644 mpmath/ctx_fp.py                     %{buildroot}%{prefix}/DSR/mpmath/ctx_fp.py
install -m 644 mpmath/ctx_iv.py                     %{buildroot}%{prefix}/DSR/mpmath/ctx_iv.py
install -m 644 mpmath/ctx_mp.py                     %{buildroot}%{prefix}/DSR/mpmath/ctx_mp.py
install -m 644 mpmath/ctx_mp_python.py              %{buildroot}%{prefix}/DSR/mpmath/ctx_mp_python.py
install -m 644 mpmath/function_docs.py              %{buildroot}%{prefix}/DSR/mpmath/function_docs.py
install -m 644 mpmath/identification.py             %{buildroot}%{prefix}/DSR/mpmath/identification.py
install -m 644 mpmath/math2.py                      %{buildroot}%{prefix}/DSR/mpmath/math2.py
install -m 644 mpmath/rational.py                   %{buildroot}%{prefix}/DSR/mpmath/rational.py
install -m 644 mpmath/usertools.py                  %{buildroot}%{prefix}/DSR/mpmath/usertools.py
install -m 644 mpmath/visualization.py              %{buildroot}%{prefix}/DSR/mpmath/visualization.py
install -m 644 mpmath/__init__.py                   %{buildroot}%{prefix}/DSR/mpmath/__init__.py
install -m 644 mpmath/calculus/approximation.py     %{buildroot}%{prefix}/DSR/mpmath/calculus/approximation.py
install -m 644 mpmath/calculus/calculus.py          %{buildroot}%{prefix}/DSR/mpmath/calculus/calculus.py
install -m 644 mpmath/calculus/differentiation.py   %{buildroot}%{prefix}/DSR/mpmath/calculus/differentiation.py
install -m 644 mpmath/calculus/extrapolation.py     %{buildroot}%{prefix}/DSR/mpmath/calculus/extrapolation.py
install -m 644 mpmath/calculus/odes.py              %{buildroot}%{prefix}/DSR/mpmath/calculus/odes.py
install -m 644 mpmath/calculus/optimization.py      %{buildroot}%{prefix}/DSR/mpmath/calculus/optimization.py
install -m 644 mpmath/calculus/polynomials.py       %{buildroot}%{prefix}/DSR/mpmath/calculus/polynomials.py
install -m 644 mpmath/calculus/quadrature.py        %{buildroot}%{prefix}/DSR/mpmath/calculus/quadrature.py
install -m 644 mpmath/calculus/__init__.py          %{buildroot}%{prefix}/DSR/mpmath/calculus/__init__.py
install -m 644 mpmath/functions/bessel.py           %{buildroot}%{prefix}/DSR/mpmath/functions/bessel.py
install -m 644 mpmath/functions/elliptic.py         %{buildroot}%{prefix}/DSR/mpmath/functions/elliptic.py
install -m 644 mpmath/functions/expintegrals.py     %{buildroot}%{prefix}/DSR/mpmath/functions/expintegrals.py
install -m 644 mpmath/functions/factorials.py       %{buildroot}%{prefix}/DSR/mpmath/functions/factorials.py
install -m 644 mpmath/functions/functions.py        %{buildroot}%{prefix}/DSR/mpmath/functions/functions.py
install -m 644 mpmath/functions/hypergeometric.py   %{buildroot}%{prefix}/DSR/mpmath/functions/hypergeometric.py
install -m 644 mpmath/functions/orthogonal.py       %{buildroot}%{prefix}/DSR/mpmath/functions/orthogonal.py
install -m 644 mpmath/functions/qfunctions.py       %{buildroot}%{prefix}/DSR/mpmath/functions/qfunctions.py
install -m 644 mpmath/functions/rszeta.py           %{buildroot}%{prefix}/DSR/mpmath/functions/rszeta.py
install -m 644 mpmath/functions/theta.py            %{buildroot}%{prefix}/DSR/mpmath/functions/theta.py
install -m 644 mpmath/functions/zeta.py             %{buildroot}%{prefix}/DSR/mpmath/functions/zeta.py
install -m 644 mpmath/functions/zetazeros.py        %{buildroot}%{prefix}/DSR/mpmath/functions/zetazeros.py
install -m 644 mpmath/functions/__init__.py         %{buildroot}%{prefix}/DSR/mpmath/functions/__init__.py
install -m 644 mpmath/libmp/backend.py              %{buildroot}%{prefix}/DSR/mpmath/libmp/backend.py
install -m 644 mpmath/libmp/gammazeta.py            %{buildroot}%{prefix}/DSR/mpmath/libmp/gammazeta.py
install -m 644 mpmath/libmp/libelefun.py            %{buildroot}%{prefix}/DSR/mpmath/libmp/libelefun.py
install -m 644 mpmath/libmp/libhyper.py             %{buildroot}%{prefix}/DSR/mpmath/libmp/libhyper.py
install -m 644 mpmath/libmp/libintmath.py           %{buildroot}%{prefix}/DSR/mpmath/libmp/libintmath.py
install -m 644 mpmath/libmp/libmpc.py               %{buildroot}%{prefix}/DSR/mpmath/libmp/libmpc.py
install -m 644 mpmath/libmp/libmpf.py               %{buildroot}%{prefix}/DSR/mpmath/libmp/libmpf.py
install -m 644 mpmath/libmp/libmpi.py               %{buildroot}%{prefix}/DSR/mpmath/libmp/libmpi.py
install -m 644 mpmath/libmp/six.py                  %{buildroot}%{prefix}/DSR/mpmath/libmp/six.py
install -m 644 mpmath/libmp/__init__.py             %{buildroot}%{prefix}/DSR/mpmath/libmp/__init__.py
install -m 644 mpmath/matrices/calculus.py          %{buildroot}%{prefix}/DSR/mpmath/matrices/calculus.py
install -m 644 mpmath/matrices/eigen.py             %{buildroot}%{prefix}/DSR/mpmath/matrices/eigen.py
install -m 644 mpmath/matrices/eigen_symmetric.py   %{buildroot}%{prefix}/DSR/mpmath/matrices/eigen_symmetric.py
install -m 644 mpmath/matrices/linalg.py            %{buildroot}%{prefix}/DSR/mpmath/matrices/linalg.py
install -m 644 mpmath/matrices/matrices.py          %{buildroot}%{prefix}/DSR/mpmath/matrices/matrices.py
install -m 644 mpmath/matrices/__init__.py          %{buildroot}%{prefix}/DSR/mpmath/matrices/__init__.py

dos2unix -q %{buildroot}%{prefix}/*



%files
%doc %{prefix}/DSR/manuals/DSR-manual.pdf
%config /etc/profile.d/dsr.sh
%{prefix}/DSR/README
%{prefix}/DSR/afix.py
%{prefix}/DSR/atomhandling.py
%{prefix}/DSR/atoms.py
%{prefix}/DSR/constants.py
%{prefix}/DSR/dbfile.py
%{prefix}/DSR/dsr.py
%{prefix}/DSR/dsr
%{prefix}/DSR/dsrparse.py
%{prefix}/DSR/export.py
%{prefix}/DSR/misc.py
%{prefix}/DSR/terminalsize.py
%{prefix}/DSR/options.py
%{prefix}/DSR/resfile.py
%{prefix}/DSR/elements.py
%{prefix}/DSR/restraints.py
%{prefix}/DSR/resi.py
%{prefix}/DSR/refine.py
%{prefix}/DSR/cf3fit.py
%{prefix}/DSR/selfupdate.py
%{prefix}/DSR/pyperclip.py
%{prefix}/DSR/dsr_db.txt
#%{prefix}/DSR/manuals/DSR-manual.pdf
%{prefix}/DSR/fit/__init__.py
%{prefix}/DSR/fit/quatfit.py
%{prefix}/DSR/fit/LICENSE.txt
%{prefix}/DSR/example/p21c.hkl
%{prefix}/DSR/example/p21c.res
%{prefix}/DSR/example/p21c_step0.res
%{prefix}/DSR/example/p21c_step1.res
%{prefix}/DSR/example/p21c_step2.ins
%{prefix}/DSR/example/p21c_step3.res
%{prefix}/DSR/example/p21c_final.res
%{prefix}/DSR/example/p21n_cf3.hkl
%{prefix}/DSR/example/p21n_cf3.res
%{prefix}/DSR/networkx/external/decorator/decorator3/_decorator3.py
%{prefix}/DSR/networkx/external/decorator/decorator3/__init__.py
%{prefix}/DSR/networkx/external/decorator/__init__.py
%{prefix}/DSR/networkx/external/decorator/decorator2/_decorator2.py
%{prefix}/DSR/networkx/external/decorator/decorator2/__init__.py
%{prefix}/DSR/networkx/external/__init__.py
%{prefix}/DSR/networkx/algorithms/richclub.py
%{prefix}/DSR/networkx/algorithms/assortativity/pairs.py
%{prefix}/DSR/networkx/algorithms/assortativity/correlation.py
%{prefix}/DSR/networkx/algorithms/assortativity/mixing.py
%{prefix}/DSR/networkx/algorithms/assortativity/__init__.py
%{prefix}/DSR/networkx/algorithms/assortativity/neighbor_degree.py
%{prefix}/DSR/networkx/algorithms/assortativity/connectivity.py
%{prefix}/DSR/networkx/algorithms/link_analysis/pagerank_alg.py
%{prefix}/DSR/networkx/algorithms/link_analysis/__init__.py
%{prefix}/DSR/networkx/algorithms/link_analysis/hits_alg.py
%{prefix}/DSR/networkx/algorithms/smetric.py
%{prefix}/DSR/networkx/algorithms/graphical.py
%{prefix}/DSR/networkx/algorithms/shortest_paths/dense.py
%{prefix}/DSR/networkx/algorithms/shortest_paths/unweighted.py
%{prefix}/DSR/networkx/algorithms/shortest_paths/weighted.py
%{prefix}/DSR/networkx/algorithms/shortest_paths/astar.py
%{prefix}/DSR/networkx/algorithms/shortest_paths/__init__.py
%{prefix}/DSR/networkx/algorithms/shortest_paths/generic.py
%{prefix}/DSR/networkx/algorithms/simple_paths.py
%{prefix}/DSR/networkx/algorithms/vitality.py
%{prefix}/DSR/networkx/algorithms/approximation/independent_set.py
%{prefix}/DSR/networkx/algorithms/approximation/dominating_set.py
%{prefix}/DSR/networkx/algorithms/approximation/matching.py
%{prefix}/DSR/networkx/algorithms/approximation/ramsey.py
%{prefix}/DSR/networkx/algorithms/approximation/vertex_cover.py
%{prefix}/DSR/networkx/algorithms/approximation/clique.py
%{prefix}/DSR/networkx/algorithms/approximation/__init__.py
%{prefix}/DSR/networkx/algorithms/chordal/chordal_alg.py
%{prefix}/DSR/networkx/algorithms/chordal/__init__.py
%{prefix}/DSR/networkx/algorithms/hierarchy.py
%{prefix}/DSR/networkx/algorithms/block.py
%{prefix}/DSR/networkx/algorithms/core.py
%{prefix}/DSR/networkx/algorithms/distance_regular.py
%{prefix}/DSR/networkx/algorithms/components/strongly_connected.py
%{prefix}/DSR/networkx/algorithms/components/weakly_connected.py
%{prefix}/DSR/networkx/algorithms/components/attracting.py
%{prefix}/DSR/networkx/algorithms/components/biconnected.py
%{prefix}/DSR/networkx/algorithms/components/__init__.py
%{prefix}/DSR/networkx/algorithms/components/connected.py
%{prefix}/DSR/networkx/algorithms/centrality/closeness.py
%{prefix}/DSR/networkx/algorithms/centrality/katz.py
%{prefix}/DSR/networkx/algorithms/centrality/betweenness.py
%{prefix}/DSR/networkx/algorithms/centrality/flow_matrix.py
%{prefix}/DSR/networkx/algorithms/centrality/load.py
%{prefix}/DSR/networkx/algorithms/centrality/current_flow_closeness.py
%{prefix}/DSR/networkx/algorithms/centrality/communicability_alg.py
%{prefix}/DSR/networkx/algorithms/centrality/eigenvector.py
%{prefix}/DSR/networkx/algorithms/centrality/degree_alg.py
%{prefix}/DSR/networkx/algorithms/centrality/current_flow_betweenness_subset.py
%{prefix}/DSR/networkx/algorithms/centrality/current_flow_betweenness.py
%{prefix}/DSR/networkx/algorithms/centrality/betweenness_subset.py
%{prefix}/DSR/networkx/algorithms/centrality/__init__.py
%{prefix}/DSR/networkx/algorithms/matching.py
%{prefix}/DSR/networkx/algorithms/dag.py
%{prefix}/DSR/networkx/algorithms/boundary.py
%{prefix}/DSR/networkx/algorithms/euler.py
%{prefix}/DSR/networkx/algorithms/flow/mincost.py
%{prefix}/DSR/networkx/algorithms/flow/maxflow.py
%{prefix}/DSR/networkx/algorithms/flow/__init__.py
%{prefix}/DSR/networkx/algorithms/distance_measures.py
%{prefix}/DSR/networkx/algorithms/isolate.py
%{prefix}/DSR/networkx/algorithms/community/kclique.py
%{prefix}/DSR/networkx/algorithms/community/__init__.py
%{prefix}/DSR/networkx/algorithms/swap.py
%{prefix}/DSR/networkx/algorithms/cycles.py
%{prefix}/DSR/networkx/algorithms/connectivity/__init__.py
%{prefix}/DSR/networkx/algorithms/connectivity/connectivity.py
%{prefix}/DSR/networkx/algorithms/connectivity/cuts.py
%{prefix}/DSR/networkx/algorithms/clique.py
%{prefix}/DSR/networkx/algorithms/__init__.py
%{prefix}/DSR/networkx/algorithms/mst.py
%{prefix}/DSR/networkx/algorithms/bipartite/centrality.py
%{prefix}/DSR/networkx/algorithms/bipartite/projection.py
%{prefix}/DSR/networkx/algorithms/bipartite/basic.py
%{prefix}/DSR/networkx/algorithms/bipartite/spectral.py
%{prefix}/DSR/networkx/algorithms/bipartite/__init__.py
%{prefix}/DSR/networkx/algorithms/bipartite/redundancy.py
%{prefix}/DSR/networkx/algorithms/bipartite/cluster.py
%{prefix}/DSR/networkx/algorithms/traversal/depth_first_search.py
%{prefix}/DSR/networkx/algorithms/traversal/__init__.py
%{prefix}/DSR/networkx/algorithms/traversal/breadth_first_search.py
%{prefix}/DSR/networkx/algorithms/mis.py
%{prefix}/DSR/networkx/algorithms/operators/unary.py
%{prefix}/DSR/networkx/algorithms/operators/all.py
%{prefix}/DSR/networkx/algorithms/operators/product.py
%{prefix}/DSR/networkx/algorithms/operators/binary.py
%{prefix}/DSR/networkx/algorithms/operators/__init__.py
%{prefix}/DSR/networkx/algorithms/isomorphism/isomorphvf2.py
%{prefix}/DSR/networkx/algorithms/isomorphism/matchhelpers.py
%{prefix}/DSR/networkx/algorithms/isomorphism/__init__.py
%{prefix}/DSR/networkx/algorithms/isomorphism/isomorph.py
%{prefix}/DSR/networkx/algorithms/isomorphism/vf2userfunc.py
%{prefix}/DSR/networkx/algorithms/cluster.py
%{prefix}/DSR/networkx/generators/classic.py
%{prefix}/DSR/networkx/generators/social.py
%{prefix}/DSR/networkx/generators/intersection.py
%{prefix}/DSR/networkx/generators/threshold.py
%{prefix}/DSR/networkx/generators/ego.py
%{prefix}/DSR/networkx/generators/atlas.py
%{prefix}/DSR/networkx/generators/hybrid.py
%{prefix}/DSR/networkx/generators/directed.py
%{prefix}/DSR/networkx/generators/random_graphs.py
%{prefix}/DSR/networkx/generators/bipartite.py
%{prefix}/DSR/networkx/generators/line.py
%{prefix}/DSR/networkx/generators/__init__.py
%{prefix}/DSR/networkx/generators/geometric.py
%{prefix}/DSR/networkx/generators/random_clustered.py
%{prefix}/DSR/networkx/generators/small.py
%{prefix}/DSR/networkx/generators/degree_seq.py
%{prefix}/DSR/networkx/generators/stochastic.py
%{prefix}/DSR/networkx/exception.py
%{prefix}/DSR/networkx/version.py
%{prefix}/DSR/networkx/release.py

%{prefix}/DSR/networkx/linalg/attrmatrix.py
%{prefix}/DSR/networkx/linalg/graphmatrix.py
%{prefix}/DSR/networkx/linalg/__init__.py
%{prefix}/DSR/networkx/linalg/laplacianmatrix.py
%{prefix}/DSR/networkx/linalg/spectrum.py
%{prefix}/DSR/networkx/readwrite/leda.py
%{prefix}/DSR/networkx/readwrite/multiline_adjlist.py
%{prefix}/DSR/networkx/readwrite/json_graph/tree.py
%{prefix}/DSR/networkx/readwrite/json_graph/serialize.py
%{prefix}/DSR/networkx/readwrite/json_graph/node_link.py
%{prefix}/DSR/networkx/readwrite/json_graph/__init__.py
%{prefix}/DSR/networkx/readwrite/json_graph/adjacency.py
%{prefix}/DSR/networkx/readwrite/edgelist.py
%{prefix}/DSR/networkx/readwrite/gml.py
%{prefix}/DSR/networkx/readwrite/pajek.py
%{prefix}/DSR/networkx/readwrite/graphml.py
%{prefix}/DSR/networkx/readwrite/gexf.py
%{prefix}/DSR/networkx/readwrite/adjlist.py
%{prefix}/DSR/networkx/readwrite/gpickle.py
%{prefix}/DSR/networkx/readwrite/__init__.py
%{prefix}/DSR/networkx/readwrite/sparsegraph6.py
%{prefix}/DSR/networkx/readwrite/nx_yaml.py
%{prefix}/DSR/networkx/readwrite/p2g.py
%{prefix}/DSR/networkx/readwrite/nx_shp.py
%{prefix}/DSR/networkx/relabel.py
%{prefix}/DSR/networkx/classes/graph.py
%{prefix}/DSR/networkx/classes/multigraph.py
%{prefix}/DSR/networkx/classes/function.py
%{prefix}/DSR/networkx/classes/digraph.py
%{prefix}/DSR/networkx/classes/__init__.py
%{prefix}/DSR/networkx/classes/multidigraph.py
%{prefix}/DSR/networkx/utils/union_find.py
%{prefix}/DSR/networkx/utils/decorators.py
%{prefix}/DSR/networkx/utils/rcm.py
%{prefix}/DSR/networkx/utils/misc.py
%{prefix}/DSR/networkx/utils/random_sequence.py
%{prefix}/DSR/networkx/utils/__init__.py
%{prefix}/DSR/networkx/drawing/nx_pydot.py
%{prefix}/DSR/networkx/drawing/__init__.py
%{prefix}/DSR/networkx/drawing/layout.py
%{prefix}/DSR/networkx/drawing/nx_pylab.py
%{prefix}/DSR/networkx/__init__.py
%{prefix}/DSR/networkx/convert.py
%{prefix}/DSR/networkx/testing/utils.py
%{prefix}/DSR/networkx/testing/__init__.py

%{prefix}/DSR/mpmath/conftest.py
%{prefix}/DSR/mpmath/ctx_base.py
%{prefix}/DSR/mpmath/ctx_fp.py
%{prefix}/DSR/mpmath/ctx_iv.py
%{prefix}/DSR/mpmath/ctx_mp.py
%{prefix}/DSR/mpmath/ctx_mp_python.py
%{prefix}/DSR/mpmath/function_docs.py
%{prefix}/DSR/mpmath/identification.py
%{prefix}/DSR/mpmath/math2.py
%{prefix}/DSR/mpmath/matrices
%{prefix}/DSR/mpmath/rational.py
%{prefix}/DSR/mpmath/usertools.py
%{prefix}/DSR/mpmath/visualization.py
%{prefix}/DSR/mpmath/__init__.py
%{prefix}/DSR/mpmath/calculus/approximation.py
%{prefix}/DSR/mpmath/calculus/calculus.py
%{prefix}/DSR/mpmath/calculus/differentiation.py
%{prefix}/DSR/mpmath/calculus/extrapolation.py
%{prefix}/DSR/mpmath/calculus/odes.py
%{prefix}/DSR/mpmath/calculus/optimization.py
%{prefix}/DSR/mpmath/calculus/polynomials.py
%{prefix}/DSR/mpmath/calculus/quadrature.py
%{prefix}/DSR/mpmath/calculus/__init__.py
%{prefix}/DSR/mpmath/functions/bessel.py
%{prefix}/DSR/mpmath/functions/elliptic.py
%{prefix}/DSR/mpmath/functions/expintegrals.py
%{prefix}/DSR/mpmath/functions/factorials.py
%{prefix}/DSR/mpmath/functions/functions.py
%{prefix}/DSR/mpmath/functions/hypergeometric.py
%{prefix}/DSR/mpmath/functions/orthogonal.py
%{prefix}/DSR/mpmath/functions/qfunctions.py
%{prefix}/DSR/mpmath/functions/rszeta.py
%{prefix}/DSR/mpmath/functions/theta.py
%{prefix}/DSR/mpmath/functions/zeta.py
%{prefix}/DSR/mpmath/functions/zetazeros.py
%{prefix}/DSR/mpmath/functions/__init__.py
%{prefix}/DSR/mpmath/libmp/backend.py
%{prefix}/DSR/mpmath/libmp/gammazeta.py
%{prefix}/DSR/mpmath/libmp/libelefun.py
%{prefix}/DSR/mpmath/libmp/libhyper.py
%{prefix}/DSR/mpmath/libmp/libintmath.py
%{prefix}/DSR/mpmath/libmp/libmpc.py
%{prefix}/DSR/mpmath/libmp/libmpf.py
%{prefix}/DSR/mpmath/libmp/libmpi.py
%{prefix}/DSR/mpmath/libmp/six.py
%{prefix}/DSR/mpmath/libmp/__init__.py
%{prefix}/DSR/mpmath/matrices/calculus.py
%{prefix}/DSR/mpmath/matrices/eigen.py
%{prefix}/DSR/mpmath/matrices/eigen_symmetric.py
%{prefix}/DSR/mpmath/matrices/linalg.py
%{prefix}/DSR/mpmath/matrices/matrices.py
%{prefix}/DSR/mpmath/matrices/__init__.py

