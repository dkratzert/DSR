Summary: DSR - A program for modelling of disordered solvents with SHELXL
Name: DSR
Provides: DSR
Packager: dkratzert@gmx.de
Version: 1.5.15
Release: 0
Requires: python, xclip
Prefix: /opt
BuildRoot: %{_tmppath}/%{name}-%{version}
BuildArch: noarch
URL: https://www.xs3.uni-freiburg.de/research/dsr
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
mkdir -p %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator3
mkdir -p %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator2
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/link_analysis/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/chordal/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/components/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/flow/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/community/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/connectivity/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/traversal/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/generators/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/linalg/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/readwrite/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/classes/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/utils/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/drawing/tests
mkdir -p %{buildroot}%{prefix}/DSR/networkx/testing/tests	
touch %{buildroot}%{prefix}/DSR/dsr_user_db.txt
mkdir -p %{buildroot}/etc/profile.d
chmod a+rw %{buildroot}%{prefix}/DSR/example

%post
echo "DSR was installed in /opt/DSR"
echo "Please restart your system to ensure that the environment variables are set."
source /etc/profile.d/dsr.sh



%install
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
install -m 666 example/p21c_step2.res %{buildroot}%{prefix}/DSR/example/p21c_step2.res
install -m 666 example/p21c_step3.res %{buildroot}%{prefix}/DSR/example/p21c_step3.res
install -m 666 example/p21c-step2.ins %{buildroot}%{prefix}/DSR/example/p21c-step2.ins

install -m 644 networkx/external/decorator/decorator3/_decorator3.py            %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator3/_decorator3.py
install -m 644 networkx/external/decorator/decorator3/__init__.py               %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator3/__init__.py
install -m 644 networkx/external/decorator/__init__.py                          %{buildroot}%{prefix}/DSR/networkx/external/decorator/__init__.py
install -m 644 networkx/external/decorator/decorator2/_decorator2.py            %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator2/_decorator2.py
install -m 644 networkx/external/decorator/decorator2/__init__.py               %{buildroot}%{prefix}/DSR/networkx/external/decorator/decorator2/__init__.py
install -m 644 networkx/external/__init__.py                                    %{buildroot}%{prefix}/DSR/networkx/external/__init__.py
install -m 644 networkx/algorithms/richclub.py                                  %{buildroot}%{prefix}/DSR/networkx/algorithms/richclub.py
install -m 644 networkx/algorithms/assortativity/tests/test_pairs.py            %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/tests/test_pairs.py
install -m 644 networkx/algorithms/assortativity/tests/test_correlation.py      %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/tests/test_correlation.py
install -m 644 networkx/algorithms/assortativity/tests/test_neighbor_degree.py  %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/tests/test_neighbor_degree.py
install -m 644 networkx/algorithms/assortativity/tests/test_mixing.py           %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/tests/test_mixing.py
install -m 644 networkx/algorithms/assortativity/tests/base_test.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/tests/base_test.py
install -m 644 networkx/algorithms/assortativity/tests/test_connectivity.py     %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/tests/test_connectivity.py
install -m 644 networkx/algorithms/assortativity/pairs.py                       %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/pairs.py
install -m 644 networkx/algorithms/assortativity/correlation.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/correlation.py
install -m 644 networkx/algorithms/assortativity/mixing.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/mixing.py
install -m 644 networkx/algorithms/assortativity/__init__.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/__init__.py
install -m 644 networkx/algorithms/assortativity/neighbor_degree.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/neighbor_degree.py                                                                                       
install -m 644 networkx/algorithms/assortativity/connectivity.py                %{buildroot}%{prefix}/DSR/networkx/algorithms/assortativity/connectivity.py                                                                                          
install -m 644 networkx/algorithms/link_analysis/tests/test_hits.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/link_analysis/tests/test_hits.py                                                                                       
install -m 644 networkx/algorithms/link_analysis/tests/test_pagerank.py         %{buildroot}%{prefix}/DSR/networkx/algorithms/link_analysis/tests/test_pagerank.py                                                                                   
install -m 644 networkx/algorithms/link_analysis/pagerank_alg.py                %{buildroot}%{prefix}/DSR/networkx/algorithms/link_analysis/pagerank_alg.py                                                                                          
install -m 644 networkx/algorithms/link_analysis/__init__.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/link_analysis/__init__.py                                                                                              
install -m 644 networkx/algorithms/link_analysis/hits_alg.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/link_analysis/hits_alg.py                                                                                              
install -m 644 networkx/algorithms/smetric.py                                   %{buildroot}%{prefix}/DSR/networkx/algorithms/smetric.py                                                                                                             
install -m 644 networkx/algorithms/graphical.py                                 %{buildroot}%{prefix}/DSR/networkx/algorithms/graphical.py                                                                                                           
install -m 644 networkx/algorithms/shortest_paths/tests/test_generic.py         %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_generic.py                                                                                   
install -m 644 networkx/algorithms/shortest_paths/tests/test_weighted.py        %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_weighted.py                                                                                  
install -m 644 networkx/algorithms/shortest_paths/tests/test_astar.py           %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_astar.py                                                                                     
install -m 644 networkx/algorithms/shortest_paths/tests/test_dense_numpy.py     %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_dense_numpy.py                                                                               
install -m 644 networkx/algorithms/shortest_paths/tests/test_dense.py           %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_dense.py                                                                                     
install -m 644 networkx/algorithms/shortest_paths/tests/test_unweighted.py      %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_unweighted.py                                                                                
install -m 644 networkx/algorithms/shortest_paths/dense.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/dense.py                                                                                                
install -m 644 networkx/algorithms/shortest_paths/unweighted.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/unweighted.py                                                                                           
install -m 644 networkx/algorithms/shortest_paths/weighted.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/weighted.py                                                                                             
install -m 644 networkx/algorithms/shortest_paths/astar.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/astar.py                                                                                                
install -m 644 networkx/algorithms/shortest_paths/__init__.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/__init__.py                                                                                             
install -m 644 networkx/algorithms/shortest_paths/generic.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/shortest_paths/generic.py                                                                                              
install -m 644 networkx/algorithms/simple_paths.py                              %{buildroot}%{prefix}/DSR/networkx/algorithms/simple_paths.py                                                                                                        
install -m 644 networkx/algorithms/vitality.py                                  %{buildroot}%{prefix}/DSR/networkx/algorithms/vitality.py                                                                                                            
install -m 644 networkx/algorithms/tests/test_mis.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_mis.py                                                                                                      
install -m 644 networkx/algorithms/tests/test_hierarchy.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_hierarchy.py                                                                                                
install -m 644 networkx/algorithms/tests/test_swap.py                           %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_swap.py                                                                                                     
install -m 644 networkx/algorithms/tests/test_simple_paths.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_simple_paths.py                                                                                             
install -m 644 networkx/algorithms/tests/test_cluster.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_cluster.py                                                                                                  
install -m 644 networkx/algorithms/tests/test_boundary.py                       %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_boundary.py                                                                                                 
install -m 644 networkx/algorithms/tests/test_distance_regular.py               %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_distance_regular.py                                                                                         
install -m 644 networkx/algorithms/tests/test_euler.py                          %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_euler.py                                                                                                    
install -m 644 networkx/algorithms/tests/test_graphical.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_graphical.py
install -m 644 networkx/algorithms/tests/test_dag.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_dag.py
install -m 644 networkx/algorithms/tests/test_smetric.py                          %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_smetric.py
install -m 644 networkx/algorithms/tests/test_richclub.py                         %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_richclub.py
install -m 644 networkx/algorithms/tests/test_distance_measures.py                %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_distance_measures.py
install -m 644 networkx/algorithms/tests/test_matching.py                         %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_matching.py
install -m 644 networkx/algorithms/tests/test_block.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_block.py
install -m 644 networkx/algorithms/tests/test_cycles.py                           %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_cycles.py
install -m 644 networkx/algorithms/tests/test_clique.py                           %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_clique.py
install -m 644 networkx/algorithms/tests/test_vitality.py                         %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_vitality.py
install -m 644 networkx/algorithms/tests/test_mst.py                              %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_mst.py
install -m 644 networkx/algorithms/tests/test_core.py                             %{buildroot}%{prefix}/DSR/networkx/algorithms/tests/test_core.py
install -m 644 networkx/algorithms/approximation/independent_set.py               %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/independent_set.py
install -m 644 networkx/algorithms/approximation/tests/test_independent_set.py    %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/tests/test_independent_set.py
install -m 644 networkx/algorithms/approximation/tests/test_vertex_cover.py       %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/tests/test_vertex_cover.py
install -m 644 networkx/algorithms/approximation/tests/test_ramsey.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/tests/test_ramsey.py
install -m 644 networkx/algorithms/approximation/tests/test_dominating_set.py     %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/tests/test_dominating_set.py
install -m 644 networkx/algorithms/approximation/tests/test_matching.py           %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/tests/test_matching.py
install -m 644 networkx/algorithms/approximation/tests/test_clique.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/tests/test_clique.py
install -m 644 networkx/algorithms/approximation/dominating_set.py                %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/dominating_set.py
install -m 644 networkx/algorithms/approximation/matching.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/matching.py
install -m 644 networkx/algorithms/approximation/ramsey.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/ramsey.py
install -m 644 networkx/algorithms/approximation/vertex_cover.py                  %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/vertex_cover.py
install -m 644 networkx/algorithms/approximation/clique.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/clique.py
install -m 644 networkx/algorithms/approximation/__init__.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/approximation/__init__.py
install -m 644 networkx/algorithms/chordal/tests/test_chordal.py                  %{buildroot}%{prefix}/DSR/networkx/algorithms/chordal/tests/test_chordal.py
install -m 644 networkx/algorithms/chordal/chordal_alg.py                         %{buildroot}%{prefix}/DSR/networkx/algorithms/chordal/chordal_alg.py
install -m 644 networkx/algorithms/chordal/__init__.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/chordal/__init__.py
install -m 644 networkx/algorithms/hierarchy.py                                   %{buildroot}%{prefix}/DSR/networkx/algorithms/hierarchy.py
install -m 644 networkx/algorithms/block.py                                       %{buildroot}%{prefix}/DSR/networkx/algorithms/block.py
install -m 644 networkx/algorithms/core.py                                        %{buildroot}%{prefix}/DSR/networkx/algorithms/core.py
install -m 644 networkx/algorithms/distance_regular.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/distance_regular.py
install -m 644 networkx/algorithms/components/tests/test_biconnected.py           %{buildroot}%{prefix}/DSR/networkx/algorithms/components/tests/test_biconnected.py
install -m 644 networkx/algorithms/components/tests/test_weakly_connected.py      %{buildroot}%{prefix}/DSR/networkx/algorithms/components/tests/test_weakly_connected.py
install -m 644 networkx/algorithms/components/tests/test_attracting.py            %{buildroot}%{prefix}/DSR/networkx/algorithms/components/tests/test_attracting.py
install -m 644 networkx/algorithms/components/tests/test_connected.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/components/tests/test_connected.py
install -m 644 networkx/algorithms/components/tests/test_strongly_connected.py    %{buildroot}%{prefix}/DSR/networkx/algorithms/components/tests/test_strongly_connected.py
install -m 644 networkx/algorithms/components/strongly_connected.py               %{buildroot}%{prefix}/DSR/networkx/algorithms/components/strongly_connected.py
install -m 644 networkx/algorithms/components/weakly_connected.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/components/weakly_connected.py
install -m 644 networkx/algorithms/components/attracting.py                       %{buildroot}%{prefix}/DSR/networkx/algorithms/components/attracting.py
install -m 644 networkx/algorithms/components/biconnected.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/components/biconnected.py
install -m 644 networkx/algorithms/components/__init__.py                         %{buildroot}%{prefix}/DSR/networkx/algorithms/components/__init__.py
install -m 644 networkx/algorithms/components/connected.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/components/connected.py
install -m 644 networkx/algorithms/centrality/closeness.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/closeness.py
install -m 644 networkx/algorithms/centrality/tests/test_load_centrality.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_load_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_betweenness_centrality_subset.py              %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_betweenness_centrality_subset.py
install -m 644 networkx/algorithms/centrality/tests/test_communicability.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_communicability.py
install -m 644 networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality_subset.py %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality_subset.py
install -m 644 networkx/algorithms/centrality/tests/test_katz_centrality.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_katz_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_current_flow_closeness.py                     %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_current_flow_closeness.py
install -m 644 networkx/algorithms/centrality/tests/test_eigenvector_centrality.py                     %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_eigenvector_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality.py        %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_closeness_centrality.py                       %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_closeness_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_betweenness_centrality.py     %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_betweenness_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_degree_centrality.py          %{buildroot}%{prefix}/DSR/networkx/algorithms/centrality/tests/test_degree_centrality.py
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
install -m 644 networkx/algorithms/flow/tests/test_mincost.py                          %{buildroot}%{prefix}/DSR/networkx/algorithms/flow/tests/test_mincost.py
install -m 644 networkx/algorithms/flow/tests/test_maxflow.py                          %{buildroot}%{prefix}/DSR/networkx/algorithms/flow/tests/test_maxflow.py
install -m 644 networkx/algorithms/flow/tests/test_maxflow_large_graph.py              %{buildroot}%{prefix}/DSR/networkx/algorithms/flow/tests/test_maxflow_large_graph.py
install -m 644 networkx/algorithms/flow/mincost.py                                     %{buildroot}%{prefix}/DSR/networkx/algorithms/flow/mincost.py
install -m 644 networkx/algorithms/flow/maxflow.py                                     %{buildroot}%{prefix}/DSR/networkx/algorithms/flow/maxflow.py
install -m 644 networkx/algorithms/flow/__init__.py                                    %{buildroot}%{prefix}/DSR/networkx/algorithms/flow/__init__.py
install -m 644 networkx/algorithms/distance_measures.py                                %{buildroot}%{prefix}/DSR/networkx/algorithms/distance_measures.py
install -m 644 networkx/algorithms/isolate.py                                          %{buildroot}%{prefix}/DSR/networkx/algorithms/isolate.py
install -m 644 networkx/algorithms/community/tests/test_kclique.py                     %{buildroot}%{prefix}/DSR/networkx/algorithms/community/tests/test_kclique.py
install -m 644 networkx/algorithms/community/kclique.py                                %{buildroot}%{prefix}/DSR/networkx/algorithms/community/kclique.py
install -m 644 networkx/algorithms/community/__init__.py                               %{buildroot}%{prefix}/DSR/networkx/algorithms/community/__init__.py
install -m 644 networkx/algorithms/swap.py                                             %{buildroot}%{prefix}/DSR/networkx/algorithms/swap.py
install -m 644 networkx/algorithms/cycles.py                                           %{buildroot}%{prefix}/DSR/networkx/algorithms/cycles.py
install -m 644 networkx/algorithms/connectivity/tests/test_cuts.py                     %{buildroot}%{prefix}/DSR/networkx/algorithms/connectivity/tests/test_cuts.py
install -m 644 networkx/algorithms/connectivity/tests/test_connectivity.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/connectivity/tests/test_connectivity.py
install -m 644 networkx/algorithms/connectivity/__init__.py                            %{buildroot}%{prefix}/DSR/networkx/algorithms/connectivity/__init__.py
install -m 644 networkx/algorithms/connectivity/connectivity.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/connectivity/connectivity.py
install -m 644 networkx/algorithms/connectivity/cuts.py                                %{buildroot}%{prefix}/DSR/networkx/algorithms/connectivity/cuts.py
install -m 644 networkx/algorithms/clique.py                                           %{buildroot}%{prefix}/DSR/networkx/algorithms/clique.py
install -m 644 networkx/algorithms/__init__.py                                         %{buildroot}%{prefix}/DSR/networkx/algorithms/__init__.py
install -m 644 networkx/algorithms/mst.py                                              %{buildroot}%{prefix}/DSR/networkx/algorithms/mst.py
install -m 644 networkx/algorithms/bipartite/centrality.py                             %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/centrality.py
install -m 644 networkx/algorithms/bipartite/projection.py                             %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/projection.py
install -m 644 networkx/algorithms/bipartite/tests/test_basic.py                       %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/tests/test_basic.py
install -m 644 networkx/algorithms/bipartite/tests/test_centrality.py                  %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/tests/test_centrality.py
install -m 644 networkx/algorithms/bipartite/tests/test_cluster.py                     %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/tests/test_cluster.py
install -m 644 networkx/algorithms/bipartite/tests/test_spectral_bipartivity.py        %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/tests/test_spectral_bipartivity.py
install -m 644 networkx/algorithms/bipartite/tests/test_project.py                     %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/tests/test_project.py
install -m 644 networkx/algorithms/bipartite/basic.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/basic.py
install -m 644 networkx/algorithms/bipartite/spectral.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/spectral.py
install -m 644 networkx/algorithms/bipartite/__init__.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/__init__.py
install -m 644 networkx/algorithms/bipartite/redundancy.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/redundancy.py
install -m 644 networkx/algorithms/bipartite/cluster.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/bipartite/cluster.py
install -m 644 networkx/algorithms/traversal/tests/test_bfs.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/traversal/tests/test_bfs.py
install -m 644 networkx/algorithms/traversal/tests/test_dfs.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/traversal/tests/test_dfs.py
install -m 644 networkx/algorithms/traversal/depth_first_search.py         %{buildroot}%{prefix}/DSR/networkx/algorithms/traversal/depth_first_search.py
install -m 644 networkx/algorithms/traversal/__init__.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/traversal/__init__.py
install -m 644 networkx/algorithms/traversal/breadth_first_search.py       %{buildroot}%{prefix}/DSR/networkx/algorithms/traversal/breadth_first_search.py
install -m 644 networkx/algorithms/mis.py                                  %{buildroot}%{prefix}/DSR/networkx/algorithms/mis.py
install -m 644 networkx/algorithms/operators/unary.py                      %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/unary.py
install -m 644 networkx/algorithms/operators/tests/test_product.py         %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/tests/test_product.py
install -m 644 networkx/algorithms/operators/tests/test_binary.py          %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/tests/test_binary.py
install -m 644 networkx/algorithms/operators/tests/test_all.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/tests/test_all.py
install -m 644 networkx/algorithms/operators/tests/test_unary.py           %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/tests/test_unary.py
install -m 644 networkx/algorithms/operators/all.py                        %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/all.py
install -m 644 networkx/algorithms/operators/product.py                    %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/product.py
install -m 644 networkx/algorithms/operators/binary.py                     %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/binary.py
install -m 644 networkx/algorithms/operators/__init__.py                   %{buildroot}%{prefix}/DSR/networkx/algorithms/operators/__init__.py
install -m 644 networkx/algorithms/isomorphism/isomorphvf2.py              %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/isomorphvf2.py
install -m 644 networkx/algorithms/isomorphism/tests/iso_r01_s80.A99       %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/tests/iso_r01_s80.A99
install -m 644 networkx/algorithms/isomorphism/tests/test_isomorphvf2.py   %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/tests/test_isomorphvf2.py
install -m 644 networkx/algorithms/isomorphism/tests/iso_r01_s80.B99       %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/tests/iso_r01_s80.B99
install -m 644 networkx/algorithms/isomorphism/tests/si2_b06_m200.B99      %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/tests/si2_b06_m200.B99
install -m 644 networkx/algorithms/isomorphism/tests/test_vf2userfunc.py   %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/tests/test_vf2userfunc.py
install -m 644 networkx/algorithms/isomorphism/tests/si2_b06_m200.A99      %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/tests/si2_b06_m200.A99
install -m 644 networkx/algorithms/isomorphism/tests/test_isomorphism.py   %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/tests/test_isomorphism.py
install -m 644 networkx/algorithms/isomorphism/matchhelpers.py             %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/matchhelpers.py
install -m 644 networkx/algorithms/isomorphism/__init__.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/__init__.py
install -m 644 networkx/algorithms/isomorphism/isomorph.py                 %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/isomorph.py
install -m 644 networkx/algorithms/isomorphism/vf2userfunc.py              %{buildroot}%{prefix}/DSR/networkx/algorithms/isomorphism/vf2userfunc.py
install -m 644 networkx/algorithms/cluster.py                              %{buildroot}%{prefix}/DSR/networkx/algorithms/cluster.py
install -m 644 networkx/tests/test_convert.py                              %{buildroot}%{prefix}/DSR/networkx/tests/test_convert.py
install -m 644 networkx/tests/test_convert_numpy.py                        %{buildroot}%{prefix}/DSR/networkx/tests/test_convert_numpy.py
install -m 644 networkx/tests/test_exceptions.py                           %{buildroot}%{prefix}/DSR/networkx/tests/test_exceptions.py
install -m 644 networkx/tests/test_relabel.py                              %{buildroot}%{prefix}/DSR/networkx/tests/test_relabel.py
install -m 644 networkx/tests/test.py                                      %{buildroot}%{prefix}/DSR/networkx/tests/test.py
install -m 644 networkx/tests/__init__.py                                  %{buildroot}%{prefix}/DSR/networkx/tests/__init__.py
install -m 644 networkx/tests/test_convert_scipy.py                        %{buildroot}%{prefix}/DSR/networkx/tests/test_convert_scipy.py
install -m 644 networkx/tests/benchmark.py                                 %{buildroot}%{prefix}/DSR/networkx/tests/benchmark.py
install -m 644 networkx/generators/tests/test_random_graphs.py             %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_random_graphs.py
install -m 644 networkx/generators/tests/test_small.py                     %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_small.py
install -m 644 networkx/generators/tests/test_geometric.py                 %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_geometric.py
install -m 644 networkx/generators/tests/test_stochastic.py                %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_stochastic.py
install -m 644 networkx/generators/tests/test_classic.py                   %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_classic.py
install -m 644 networkx/generators/tests/test_threshold.py                 %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_threshold.py
install -m 644 networkx/generators/tests/test_atlas.py                     %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_atlas.py
install -m 644 networkx/generators/tests/test_line.py                      %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_line.py
install -m 644 networkx/generators/tests/test_random_clustered.py          %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_random_clustered.py
install -m 644 networkx/generators/tests/test_hybrid.py                    %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_hybrid.py
install -m 644 networkx/generators/tests/test_directed.py                  %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_directed.py
install -m 644 networkx/generators/tests/test_degree_seq.py                %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_degree_seq.py
install -m 644 networkx/generators/tests/test_bipartite.py                 %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_bipartite.py
install -m 644 networkx/generators/tests/test_intersection.py              %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_intersection.py
install -m 644 networkx/generators/tests/test_ego.py                       %{buildroot}%{prefix}/DSR/networkx/generators/tests/test_ego.py
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
install -m 644 networkx/linalg/tests/test_spectrum.py                      %{buildroot}%{prefix}/DSR/networkx/linalg/tests/test_spectrum.py
install -m 644 networkx/linalg/tests/test_laplacian.py                     %{buildroot}%{prefix}/DSR/networkx/linalg/tests/test_laplacian.py
install -m 644 networkx/linalg/tests/test_graphmatrix.py                   %{buildroot}%{prefix}/DSR/networkx/linalg/tests/test_graphmatrix.py
install -m 644 networkx/linalg/attrmatrix.py                               %{buildroot}%{prefix}/DSR/networkx/linalg/attrmatrix.py
install -m 644 networkx/linalg/graphmatrix.py                              %{buildroot}%{prefix}/DSR/networkx/linalg/graphmatrix.py
install -m 644 networkx/linalg/__init__.py                                 %{buildroot}%{prefix}/DSR/networkx/linalg/__init__.py
install -m 644 networkx/linalg/laplacianmatrix.py                          %{buildroot}%{prefix}/DSR/networkx/linalg/laplacianmatrix.py
install -m 644 networkx/linalg/spectrum.py                                 %{buildroot}%{prefix}/DSR/networkx/linalg/spectrum.py
install -m 644 networkx/readwrite/leda.py                                  %{buildroot}%{prefix}/DSR/networkx/readwrite/leda.py
install -m 644 networkx/readwrite/tests/test_p2g.py                        %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_p2g.py
install -m 644 networkx/readwrite/tests/test_leda.py                       %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_leda.py
install -m 644 networkx/readwrite/tests/test_shp.py                        %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_shp.py
install -m 644 networkx/readwrite/tests/test_gpickle.py                    %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_gpickle.py
install -m 644 networkx/readwrite/tests/test_yaml.py                       %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_yaml.py
install -m 644 networkx/readwrite/tests/test_graphml.py                    %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_graphml.py
install -m 644 networkx/readwrite/tests/test_gexf.py                       %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_gexf.py
install -m 644 networkx/readwrite/tests/test_sparsegraph6.py               %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_sparsegraph6.py
install -m 644 networkx/readwrite/tests/test_gml.py                        %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_gml.py
install -m 644 networkx/readwrite/tests/test_adjlist.py                    %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_adjlist.py
install -m 644 networkx/readwrite/tests/test_edgelist.py                   %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_edgelist.py
install -m 644 networkx/readwrite/tests/test_pajek.py                      %{buildroot}%{prefix}/DSR/networkx/readwrite/tests/test_pajek.py
install -m 644 networkx/readwrite/multiline_adjlist.py                     %{buildroot}%{prefix}/DSR/networkx/readwrite/multiline_adjlist.py
install -m 644 networkx/readwrite/json_graph/tests/test_tree.py            %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph/tests/test_tree.py
install -m 644 networkx/readwrite/json_graph/tests/test_adjacency.py       %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph/tests/test_adjacency.py
install -m 644 networkx/readwrite/json_graph/tests/test_serialize.py       %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph/tests/test_serialize.py
install -m 644 networkx/readwrite/json_graph/tests/test_node_link.py       %{buildroot}%{prefix}/DSR/networkx/readwrite/json_graph/tests/test_node_link.py
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
install -m 644 networkx/classes/tests/test_multigraph.py                   %{buildroot}%{prefix}/DSR/networkx/classes/tests/test_multigraph.py
install -m 644 networkx/classes/tests/test_graph.py                        %{buildroot}%{prefix}/DSR/networkx/classes/tests/test_graph.py
install -m 644 networkx/classes/tests/test_multidigraph.py                 %{buildroot}%{prefix}/DSR/networkx/classes/tests/test_multidigraph.py
install -m 644 networkx/classes/tests/test_digraph.py                      %{buildroot}%{prefix}/DSR/networkx/classes/tests/test_digraph.py
install -m 644 networkx/classes/tests/test_graph_historical.py             %{buildroot}%{prefix}/DSR/networkx/classes/tests/test_graph_historical.py
install -m 644 networkx/classes/tests/historical_tests.py                  %{buildroot}%{prefix}/DSR/networkx/classes/tests/historical_tests.py
install -m 644 networkx/classes/tests/test_digraph_historical.py           %{buildroot}%{prefix}/DSR/networkx/classes/tests/test_digraph_historical.py
install -m 644 networkx/classes/tests/test_function.py                     %{buildroot}%{prefix}/DSR/networkx/classes/tests/test_function.py
install -m 644 networkx/classes/graph.py                                   %{buildroot}%{prefix}/DSR/networkx/classes/graph.py
install -m 644 networkx/classes/multigraph.py                              %{buildroot}%{prefix}/DSR/networkx/classes/multigraph.py
install -m 644 networkx/classes/function.py                                %{buildroot}%{prefix}/DSR/networkx/classes/function.py
install -m 644 networkx/classes/digraph.py                                 %{buildroot}%{prefix}/DSR/networkx/classes/digraph.py
install -m 644 networkx/classes/__init__.py                                %{buildroot}%{prefix}/DSR/networkx/classes/__init__.py
install -m 644 networkx/classes/multidigraph.py                            %{buildroot}%{prefix}/DSR/networkx/classes/multidigraph.py
install -m 644 networkx/utils/tests/test_misc.py                           %{buildroot}%{prefix}/DSR/networkx/utils/tests/test_misc.py
install -m 644 networkx/utils/tests/test_decorators.py                     %{buildroot}%{prefix}/DSR/networkx/utils/tests/test_decorators.py
install -m 644 networkx/utils/tests/test_random_sequence.py                %{buildroot}%{prefix}/DSR/networkx/utils/tests/test_random_sequence.py
install -m 644 networkx/utils/tests/test_rcm.py                            %{buildroot}%{prefix}/DSR/networkx/utils/tests/test_rcm.py
install -m 644 networkx/utils/tests/test.txt                               %{buildroot}%{prefix}/DSR/networkx/utils/tests/test.txt
install -m 644 networkx/utils/union_find.py                                %{buildroot}%{prefix}/DSR/networkx/utils/union_find.py
install -m 644 networkx/utils/decorators.py                                %{buildroot}%{prefix}/DSR/networkx/utils/decorators.py
install -m 644 networkx/utils/rcm.py                                       %{buildroot}%{prefix}/DSR/networkx/utils/rcm.py
install -m 644 networkx/utils/misc.py                                      %{buildroot}%{prefix}/DSR/networkx/utils/misc.py
install -m 644 networkx/utils/random_sequence.py                           %{buildroot}%{prefix}/DSR/networkx/utils/random_sequence.py
install -m 644 networkx/utils/__init__.py                                  %{buildroot}%{prefix}/DSR/networkx/utils/__init__.py
install -m 644 networkx/drawing/tests/test_agraph.py                       %{buildroot}%{prefix}/DSR/networkx/drawing/tests/test_agraph.py
install -m 644 networkx/drawing/tests/test_pydot.py                        %{buildroot}%{prefix}/DSR/networkx/drawing/tests/test_pydot.py
install -m 644 networkx/drawing/tests/test_layout.py                       %{buildroot}%{prefix}/DSR/networkx/drawing/tests/test_layout.py
install -m 644 networkx/drawing/tests/test_pylab.py                        %{buildroot}%{prefix}/DSR/networkx/drawing/tests/test_pylab.py
install -m 644 networkx/drawing/nx_agraph.py                               %{buildroot}%{prefix}/DSR/networkx/drawing/nx_agraph.py
install -m 644 networkx/drawing/nx_pydot.py                                %{buildroot}%{prefix}/DSR/networkx/drawing/nx_pydot.py
install -m 644 networkx/drawing/__init__.py                                %{buildroot}%{prefix}/DSR/networkx/drawing/__init__.py
install -m 644 networkx/drawing/layout.py                                  %{buildroot}%{prefix}/DSR/networkx/drawing/layout.py
install -m 644 networkx/drawing/nx_pylab.py                                %{buildroot}%{prefix}/DSR/networkx/drawing/nx_pylab.py
install -m 644 networkx/__init__.py                                        %{buildroot}%{prefix}/DSR/networkx/__init__.py
install -m 644 networkx/convert.py                                         %{buildroot}%{prefix}/DSR/networkx/convert.py
install -m 644 networkx/testing/tests/test_utils.py                        %{buildroot}%{prefix}/DSR/networkx/testing/tests/test_utils.py
install -m 644 networkx/testing/utils.py                                   %{buildroot}%{prefix}/DSR/networkx/testing/utils.py
install -m 644 networkx/testing/__init__.py                                %{buildroot}%{prefix}/DSR/networkx/testing/__init__.py

dos2unix -q %{buildroot}%{prefix}/*



%files
%doc %{prefix}/DSR/manuals/DSR-manual.pdf
%config /etc/profile.d/dsr.sh
#%{prefix}/DSR/dsr
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
%{prefix}/DSR/restraints.py
%{prefix}/DSR/resi.py
%{prefix}/DSR/refine.py
%{prefix}/DSR/pyperclip.py
%{prefix}/DSR/dsr_db.txt
%{prefix}/DSR/dsr_user_db.txt
#%{prefix}/DSR/manuals/DSR-manual.pdf
%{prefix}/DSR/example/p21c.hkl
%{prefix}/DSR/example/p21c.res
%{prefix}/DSR/example/p21c_step0.res
%{prefix}/DSR/example/p21c_step1.res
%{prefix}/DSR/example/p21c_step2.res
%{prefix}/DSR/example/p21c_step3.res
%{prefix}/DSR/example/p21c-step2.ins
%{prefix}/DSR/networkx/external/decorator/decorator3/_decorator3.py
%{prefix}/DSR/networkx/external/decorator/decorator3/__init__.py
%{prefix}/DSR/networkx/external/decorator/__init__.py
%{prefix}/DSR/networkx/external/decorator/decorator2/_decorator2.py
%{prefix}/DSR/networkx/external/decorator/decorator2/__init__.py
%{prefix}/DSR/networkx/external/__init__.py
%{prefix}/DSR/networkx/algorithms/richclub.py
%{prefix}/DSR/networkx/algorithms/assortativity/tests/test_pairs.py
%{prefix}/DSR/networkx/algorithms/assortativity/tests/test_correlation.py
%{prefix}/DSR/networkx/algorithms/assortativity/tests/test_neighbor_degree.py
%{prefix}/DSR/networkx/algorithms/assortativity/tests/test_mixing.py
%{prefix}/DSR/networkx/algorithms/assortativity/tests/base_test.py
%{prefix}/DSR/networkx/algorithms/assortativity/tests/test_connectivity.py
%{prefix}/DSR/networkx/algorithms/assortativity/pairs.py
%{prefix}/DSR/networkx/algorithms/assortativity/correlation.py
%{prefix}/DSR/networkx/algorithms/assortativity/mixing.py
%{prefix}/DSR/networkx/algorithms/assortativity/__init__.py
%{prefix}/DSR/networkx/algorithms/assortativity/neighbor_degree.py                                                                                            
%{prefix}/DSR/networkx/algorithms/assortativity/connectivity.py                                                                                               
%{prefix}/DSR/networkx/algorithms/link_analysis/tests/test_hits.py                                                                                            
%{prefix}/DSR/networkx/algorithms/link_analysis/tests/test_pagerank.py                                                                                        
%{prefix}/DSR/networkx/algorithms/link_analysis/pagerank_alg.py                                                                                               
%{prefix}/DSR/networkx/algorithms/link_analysis/__init__.py                                                                                                   
%{prefix}/DSR/networkx/algorithms/link_analysis/hits_alg.py                                                                                                   
%{prefix}/DSR/networkx/algorithms/smetric.py                                                                                                                  
%{prefix}/DSR/networkx/algorithms/graphical.py                                                                                                                
%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_generic.py                                                                                        
%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_weighted.py                                                                                       
%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_astar.py                                                                                          
%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_dense_numpy.py                                                                                    
%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_dense.py                                                                                          
%{prefix}/DSR/networkx/algorithms/shortest_paths/tests/test_unweighted.py                                                                                     
%{prefix}/DSR/networkx/algorithms/shortest_paths/dense.py                                                                                                     
%{prefix}/DSR/networkx/algorithms/shortest_paths/unweighted.py                                                                                                
%{prefix}/DSR/networkx/algorithms/shortest_paths/weighted.py                                                                                                  
%{prefix}/DSR/networkx/algorithms/shortest_paths/astar.py                                                                                                     
%{prefix}/DSR/networkx/algorithms/shortest_paths/__init__.py                                                                                                  
%{prefix}/DSR/networkx/algorithms/shortest_paths/generic.py                                                                                                   
%{prefix}/DSR/networkx/algorithms/simple_paths.py                                                                                                             
%{prefix}/DSR/networkx/algorithms/vitality.py                                                                                                                 
%{prefix}/DSR/networkx/algorithms/tests/test_mis.py                                                                                                           
%{prefix}/DSR/networkx/algorithms/tests/test_hierarchy.py                                                                                                     
%{prefix}/DSR/networkx/algorithms/tests/test_swap.py                                                                                                          
%{prefix}/DSR/networkx/algorithms/tests/test_simple_paths.py                                                                                                  
%{prefix}/DSR/networkx/algorithms/tests/test_cluster.py                                                                                                       
%{prefix}/DSR/networkx/algorithms/tests/test_boundary.py                                                                                                      
%{prefix}/DSR/networkx/algorithms/tests/test_distance_regular.py                                                                                              
%{prefix}/DSR/networkx/algorithms/tests/test_euler.py                                                                                                         
%{prefix}/DSR/networkx/algorithms/tests/test_graphical.py
%{prefix}/DSR/networkx/algorithms/tests/test_dag.py
%{prefix}/DSR/networkx/algorithms/tests/test_smetric.py
%{prefix}/DSR/networkx/algorithms/tests/test_richclub.py
%{prefix}/DSR/networkx/algorithms/tests/test_distance_measures.py
%{prefix}/DSR/networkx/algorithms/tests/test_matching.py
%{prefix}/DSR/networkx/algorithms/tests/test_block.py
%{prefix}/DSR/networkx/algorithms/tests/test_cycles.py
%{prefix}/DSR/networkx/algorithms/tests/test_clique.py
%{prefix}/DSR/networkx/algorithms/tests/test_vitality.py
%{prefix}/DSR/networkx/algorithms/tests/test_mst.py
%{prefix}/DSR/networkx/algorithms/tests/test_core.py
%{prefix}/DSR/networkx/algorithms/approximation/independent_set.py
%{prefix}/DSR/networkx/algorithms/approximation/tests/test_independent_set.py
%{prefix}/DSR/networkx/algorithms/approximation/tests/test_vertex_cover.py
%{prefix}/DSR/networkx/algorithms/approximation/tests/test_ramsey.py
%{prefix}/DSR/networkx/algorithms/approximation/tests/test_dominating_set.py
%{prefix}/DSR/networkx/algorithms/approximation/tests/test_matching.py
%{prefix}/DSR/networkx/algorithms/approximation/tests/test_clique.py
%{prefix}/DSR/networkx/algorithms/approximation/dominating_set.py
%{prefix}/DSR/networkx/algorithms/approximation/matching.py
%{prefix}/DSR/networkx/algorithms/approximation/ramsey.py
%{prefix}/DSR/networkx/algorithms/approximation/vertex_cover.py
%{prefix}/DSR/networkx/algorithms/approximation/clique.py
%{prefix}/DSR/networkx/algorithms/approximation/__init__.py
%{prefix}/DSR/networkx/algorithms/chordal/tests/test_chordal.py
%{prefix}/DSR/networkx/algorithms/chordal/chordal_alg.py
%{prefix}/DSR/networkx/algorithms/chordal/__init__.py
%{prefix}/DSR/networkx/algorithms/hierarchy.py
%{prefix}/DSR/networkx/algorithms/block.py
%{prefix}/DSR/networkx/algorithms/core.py
%{prefix}/DSR/networkx/algorithms/distance_regular.py
%{prefix}/DSR/networkx/algorithms/components/tests/test_biconnected.py
%{prefix}/DSR/networkx/algorithms/components/tests/test_weakly_connected.py
%{prefix}/DSR/networkx/algorithms/components/tests/test_attracting.py
%{prefix}/DSR/networkx/algorithms/components/tests/test_connected.py
%{prefix}/DSR/networkx/algorithms/components/tests/test_strongly_connected.py
%{prefix}/DSR/networkx/algorithms/components/strongly_connected.py
%{prefix}/DSR/networkx/algorithms/components/weakly_connected.py
%{prefix}/DSR/networkx/algorithms/components/attracting.py
%{prefix}/DSR/networkx/algorithms/components/biconnected.py
%{prefix}/DSR/networkx/algorithms/components/__init__.py
%{prefix}/DSR/networkx/algorithms/components/connected.py
%{prefix}/DSR/networkx/algorithms/centrality/closeness.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_load_centrality.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_betweenness_centrality_subset.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_communicability.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality_subset.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_katz_centrality.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_current_flow_closeness.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_eigenvector_centrality.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_closeness_centrality.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_betweenness_centrality.py
%{prefix}/DSR/networkx/algorithms/centrality/tests/test_degree_centrality.py
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
%{prefix}/DSR/networkx/algorithms/flow/tests/test_mincost.py
%{prefix}/DSR/networkx/algorithms/flow/tests/test_maxflow.py
%{prefix}/DSR/networkx/algorithms/flow/tests/test_maxflow_large_graph.py
%{prefix}/DSR/networkx/algorithms/flow/mincost.py
%{prefix}/DSR/networkx/algorithms/flow/maxflow.py
%{prefix}/DSR/networkx/algorithms/flow/__init__.py
%{prefix}/DSR/networkx/algorithms/distance_measures.py
%{prefix}/DSR/networkx/algorithms/isolate.py
%{prefix}/DSR/networkx/algorithms/community/tests/test_kclique.py
%{prefix}/DSR/networkx/algorithms/community/kclique.py
%{prefix}/DSR/networkx/algorithms/community/__init__.py
%{prefix}/DSR/networkx/algorithms/swap.py
%{prefix}/DSR/networkx/algorithms/cycles.py
%{prefix}/DSR/networkx/algorithms/connectivity/tests/test_cuts.py
%{prefix}/DSR/networkx/algorithms/connectivity/tests/test_connectivity.py
%{prefix}/DSR/networkx/algorithms/connectivity/__init__.py
%{prefix}/DSR/networkx/algorithms/connectivity/connectivity.py
%{prefix}/DSR/networkx/algorithms/connectivity/cuts.py
%{prefix}/DSR/networkx/algorithms/clique.py
%{prefix}/DSR/networkx/algorithms/__init__.py
%{prefix}/DSR/networkx/algorithms/mst.py
%{prefix}/DSR/networkx/algorithms/bipartite/centrality.py
%{prefix}/DSR/networkx/algorithms/bipartite/projection.py
%{prefix}/DSR/networkx/algorithms/bipartite/tests/test_basic.py
%{prefix}/DSR/networkx/algorithms/bipartite/tests/test_centrality.py
%{prefix}/DSR/networkx/algorithms/bipartite/tests/test_cluster.py
%{prefix}/DSR/networkx/algorithms/bipartite/tests/test_spectral_bipartivity.py
%{prefix}/DSR/networkx/algorithms/bipartite/tests/test_project.py
%{prefix}/DSR/networkx/algorithms/bipartite/basic.py
%{prefix}/DSR/networkx/algorithms/bipartite/spectral.py
%{prefix}/DSR/networkx/algorithms/bipartite/__init__.py
%{prefix}/DSR/networkx/algorithms/bipartite/redundancy.py
%{prefix}/DSR/networkx/algorithms/bipartite/cluster.py
%{prefix}/DSR/networkx/algorithms/traversal/tests/test_bfs.py
%{prefix}/DSR/networkx/algorithms/traversal/tests/test_dfs.py
%{prefix}/DSR/networkx/algorithms/traversal/depth_first_search.py
%{prefix}/DSR/networkx/algorithms/traversal/__init__.py
%{prefix}/DSR/networkx/algorithms/traversal/breadth_first_search.py
%{prefix}/DSR/networkx/algorithms/mis.py
%{prefix}/DSR/networkx/algorithms/operators/unary.py
%{prefix}/DSR/networkx/algorithms/operators/tests/test_product.py
%{prefix}/DSR/networkx/algorithms/operators/tests/test_binary.py
%{prefix}/DSR/networkx/algorithms/operators/tests/test_all.py
%{prefix}/DSR/networkx/algorithms/operators/tests/test_unary.py
%{prefix}/DSR/networkx/algorithms/operators/all.py
%{prefix}/DSR/networkx/algorithms/operators/product.py
%{prefix}/DSR/networkx/algorithms/operators/binary.py
%{prefix}/DSR/networkx/algorithms/operators/__init__.py
%{prefix}/DSR/networkx/algorithms/isomorphism/isomorphvf2.py
%{prefix}/DSR/networkx/algorithms/isomorphism/tests/iso_r01_s80.A99
%{prefix}/DSR/networkx/algorithms/isomorphism/tests/test_isomorphvf2.py
%{prefix}/DSR/networkx/algorithms/isomorphism/tests/iso_r01_s80.B99
%{prefix}/DSR/networkx/algorithms/isomorphism/tests/si2_b06_m200.B99
%{prefix}/DSR/networkx/algorithms/isomorphism/tests/test_vf2userfunc.py
%{prefix}/DSR/networkx/algorithms/isomorphism/tests/si2_b06_m200.A99
%{prefix}/DSR/networkx/algorithms/isomorphism/tests/test_isomorphism.py
%{prefix}/DSR/networkx/algorithms/isomorphism/matchhelpers.py
%{prefix}/DSR/networkx/algorithms/isomorphism/__init__.py
%{prefix}/DSR/networkx/algorithms/isomorphism/isomorph.py
%{prefix}/DSR/networkx/algorithms/isomorphism/vf2userfunc.py
%{prefix}/DSR/networkx/algorithms/cluster.py
%{prefix}/DSR/networkx/tests/test_convert.py
%{prefix}/DSR/networkx/tests/test_convert_numpy.py
%{prefix}/DSR/networkx/tests/test_exceptions.py
%{prefix}/DSR/networkx/tests/test_relabel.py
%{prefix}/DSR/networkx/tests/test.py
%{prefix}/DSR/networkx/tests/__init__.py
%{prefix}/DSR/networkx/tests/test_convert_scipy.py
%{prefix}/DSR/networkx/tests/benchmark.py
%{prefix}/DSR/networkx/generators/tests/test_random_graphs.py
%{prefix}/DSR/networkx/generators/tests/test_small.py
%{prefix}/DSR/networkx/generators/tests/test_geometric.py
%{prefix}/DSR/networkx/generators/tests/test_stochastic.py
%{prefix}/DSR/networkx/generators/tests/test_classic.py
%{prefix}/DSR/networkx/generators/tests/test_threshold.py
%{prefix}/DSR/networkx/generators/tests/test_atlas.py
%{prefix}/DSR/networkx/generators/tests/test_line.py
%{prefix}/DSR/networkx/generators/tests/test_random_clustered.py
%{prefix}/DSR/networkx/generators/tests/test_hybrid.py
%{prefix}/DSR/networkx/generators/tests/test_directed.py
%{prefix}/DSR/networkx/generators/tests/test_degree_seq.py
%{prefix}/DSR/networkx/generators/tests/test_bipartite.py
%{prefix}/DSR/networkx/generators/tests/test_intersection.py
%{prefix}/DSR/networkx/generators/tests/test_ego.py
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
%{prefix}/DSR/networkx/linalg/tests/test_spectrum.py
%{prefix}/DSR/networkx/linalg/tests/test_laplacian.py
%{prefix}/DSR/networkx/linalg/tests/test_graphmatrix.py
%{prefix}/DSR/networkx/linalg/attrmatrix.py
%{prefix}/DSR/networkx/linalg/graphmatrix.py
%{prefix}/DSR/networkx/linalg/__init__.py
%{prefix}/DSR/networkx/linalg/laplacianmatrix.py
%{prefix}/DSR/networkx/linalg/spectrum.py
%{prefix}/DSR/networkx/readwrite/leda.py
%{prefix}/DSR/networkx/readwrite/tests/test_p2g.py
%{prefix}/DSR/networkx/readwrite/tests/test_leda.py
%{prefix}/DSR/networkx/readwrite/tests/test_shp.py
%{prefix}/DSR/networkx/readwrite/tests/test_gpickle.py
%{prefix}/DSR/networkx/readwrite/tests/test_yaml.py
%{prefix}/DSR/networkx/readwrite/tests/test_graphml.py
%{prefix}/DSR/networkx/readwrite/tests/test_gexf.py
%{prefix}/DSR/networkx/readwrite/tests/test_sparsegraph6.py
%{prefix}/DSR/networkx/readwrite/tests/test_gml.py
%{prefix}/DSR/networkx/readwrite/tests/test_adjlist.py
%{prefix}/DSR/networkx/readwrite/tests/test_edgelist.py
%{prefix}/DSR/networkx/readwrite/tests/test_pajek.py
%{prefix}/DSR/networkx/readwrite/multiline_adjlist.py
%{prefix}/DSR/networkx/readwrite/json_graph/tests/test_tree.py
%{prefix}/DSR/networkx/readwrite/json_graph/tests/test_adjacency.py
%{prefix}/DSR/networkx/readwrite/json_graph/tests/test_serialize.py
%{prefix}/DSR/networkx/readwrite/json_graph/tests/test_node_link.py
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
%{prefix}/DSR/networkx/classes/tests/test_multigraph.py
%{prefix}/DSR/networkx/classes/tests/test_graph.py
%{prefix}/DSR/networkx/classes/tests/test_multidigraph.py
%{prefix}/DSR/networkx/classes/tests/test_digraph.py
%{prefix}/DSR/networkx/classes/tests/test_graph_historical.py
%{prefix}/DSR/networkx/classes/tests/historical_tests.py
%{prefix}/DSR/networkx/classes/tests/test_digraph_historical.py
%{prefix}/DSR/networkx/classes/tests/test_function.py
%{prefix}/DSR/networkx/classes/graph.py
%{prefix}/DSR/networkx/classes/multigraph.py
%{prefix}/DSR/networkx/classes/function.py
%{prefix}/DSR/networkx/classes/digraph.py
%{prefix}/DSR/networkx/classes/__init__.py
%{prefix}/DSR/networkx/classes/multidigraph.py
%{prefix}/DSR/networkx/utils/tests/test_misc.py
%{prefix}/DSR/networkx/utils/tests/test_decorators.py
%{prefix}/DSR/networkx/utils/tests/test_random_sequence.py
%{prefix}/DSR/networkx/utils/tests/test_rcm.py
%{prefix}/DSR/networkx/utils/tests/test.txt
%{prefix}/DSR/networkx/utils/union_find.py
%{prefix}/DSR/networkx/utils/decorators.py
%{prefix}/DSR/networkx/utils/rcm.py
%{prefix}/DSR/networkx/utils/misc.py
%{prefix}/DSR/networkx/utils/random_sequence.py
%{prefix}/DSR/networkx/utils/__init__.py
%{prefix}/DSR/networkx/drawing/tests/test_agraph.py
%{prefix}/DSR/networkx/drawing/tests/test_pydot.py
%{prefix}/DSR/networkx/drawing/tests/test_layout.py
%{prefix}/DSR/networkx/drawing/tests/test_pylab.py
%{prefix}/DSR/networkx/drawing/nx_agraph.py
%{prefix}/DSR/networkx/drawing/nx_pydot.py
%{prefix}/DSR/networkx/drawing/__init__.py
%{prefix}/DSR/networkx/drawing/layout.py
%{prefix}/DSR/networkx/drawing/nx_pylab.py
%{prefix}/DSR/networkx/__init__.py
%{prefix}/DSR/networkx/convert.py
%{prefix}/DSR/networkx/testing/tests/test_utils.py
%{prefix}/DSR/networkx/testing/utils.py
%{prefix}/DSR/networkx/testing/__init__.py