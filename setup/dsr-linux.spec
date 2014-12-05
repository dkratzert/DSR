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
mkdir %{buildroot}
mkdir %{buildroot}/DSR
mkdir %{buildroot}/DSR/networkx
mkdir %{buildroot}/DSR/manuals
mkdir %{buildroot}/DSR/example
mkdir %{buildroot}/DSR/setup
mkdir -p %{buildroot}/etc/profile.d
touch %{buildroot}/DSR/dsr_user_db.txt
find $RPM_BUILD_ROOT -type f -print > $RPM_BUILD_DIR/tmp-filelist

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

install -m 644 networkx/external/decorator/decorator3/_decorator3.py                                       %{buildroot}/DSR/networkx/external/decorator/decorator3/_decorator3.py
install -m 644 networkx/external/decorator/decorator3/__init__.py                                          %{buildroot}/DSR/networkx/external/decorator/decorator3/__init__.py
install -m 644 networkx/external/decorator/__init__.py                                                     %{buildroot}/DSR/networkx/external/decorator/__init__.py
install -m 644 networkx/external/decorator/decorator2/_decorator2.py                                       %{buildroot}/DSR/networkx/external/decorator/decorator2/_decorator2.py
install -m 644 networkx/external/decorator/decorator2/__init__.py                                          %{buildroot}/DSR/networkx/external/decorator/decorator2/__init__.py
install -m 644 networkx/external/__init__.py                                                               %{buildroot}/DSR/networkx/external/__init__.py
install -m 644 networkx/algorithms/richclub.py                                                             %{buildroot}/DSR/networkx/algorithms/richclub.py
install -m 644 networkx/algorithms/assortativity/tests/test_pairs.py                                       %{buildroot}/DSR/networkx/algorithms/assortativity/tests/test_pairs.py
install -m 644 networkx/algorithms/assortativity/tests/test_correlation.py                                 %{buildroot}/DSR/networkx/algorithms/assortativity/tests/test_correlation.py
install -m 644 networkx/algorithms/assortativity/tests/test_neighbor_degree.py                             %{buildroot}/DSR/networkx/algorithms/assortativity/tests/test_neighbor_degree.py
install -m 644 networkx/algorithms/assortativity/tests/test_mixing.py                                      %{buildroot}/DSR/networkx/algorithms/assortativity/tests/test_mixing.py
install -m 644 networkx/algorithms/assortativity/tests/base_test.py                                        %{buildroot}/DSR/networkx/algorithms/assortativity/tests/base_test.py
install -m 644 networkx/algorithms/assortativity/tests/test_connectivity.py                                %{buildroot}/DSR/networkx/algorithms/assortativity/tests/test_connectivity.py
install -m 644 networkx/algorithms/assortativity/pairs.py                                                  %{buildroot}/DSR/networkx/algorithms/assortativity/pairs.py
install -m 644 networkx/algorithms/assortativity/correlation.py                                            %{buildroot}/DSR/networkx/algorithms/assortativity/correlation.py
install -m 644 networkx/algorithms/assortativity/mixing.py                                                 %{buildroot}/DSR/networkx/algorithms/assortativity/mixing.py
install -m 644 networkx/algorithms/assortativity/__init__.py                                               %{buildroot}/DSR/networkx/algorithms/assortativity/__init__.py
install -m 644 networkx/algorithms/assortativity/neighbor_degree.py                                        %{buildroot}/DSR/networkx/algorithms/assortativity/neighbor_degree.py                                                                                       
install -m 644 networkx/algorithms/assortativity/connectivity.py                                           %{buildroot}/DSR/networkx/algorithms/assortativity/connectivity.py                                                                                          
install -m 644 networkx/algorithms/link_analysis/tests/test_hits.py                                        %{buildroot}/DSR/networkx/algorithms/link_analysis/tests/test_hits.py                                                                                       
install -m 644 networkx/algorithms/link_analysis/tests/test_pagerank.py                                    %{buildroot}/DSR/networkx/algorithms/link_analysis/tests/test_pagerank.py                                                                                   
install -m 644 networkx/algorithms/link_analysis/pagerank_alg.py                                           %{buildroot}/DSR/networkx/algorithms/link_analysis/pagerank_alg.py                                                                                          
install -m 644 networkx/algorithms/link_analysis/__init__.py                                               %{buildroot}/DSR/networkx/algorithms/link_analysis/__init__.py                                                                                              
install -m 644 networkx/algorithms/link_analysis/hits_alg.py                                               %{buildroot}/DSR/networkx/algorithms/link_analysis/hits_alg.py                                                                                              
install -m 644 networkx/algorithms/smetric.py                                                              %{buildroot}/DSR/networkx/algorithms/smetric.py                                                                                                             
install -m 644 networkx/algorithms/graphical.py                                                            %{buildroot}/DSR/networkx/algorithms/graphical.py                                                                                                           
install -m 644 networkx/algorithms/shortest_paths/tests/test_generic.py                                    %{buildroot}/DSR/networkx/algorithms/shortest_paths/tests/test_generic.py                                                                                   
install -m 644 networkx/algorithms/shortest_paths/tests/test_weighted.py                                   %{buildroot}/DSR/networkx/algorithms/shortest_paths/tests/test_weighted.py                                                                                  
install -m 644 networkx/algorithms/shortest_paths/tests/test_astar.py                                      %{buildroot}/DSR/networkx/algorithms/shortest_paths/tests/test_astar.py                                                                                     
install -m 644 networkx/algorithms/shortest_paths/tests/test_dense_numpy.py                                %{buildroot}/DSR/networkx/algorithms/shortest_paths/tests/test_dense_numpy.py                                                                               
install -m 644 networkx/algorithms/shortest_paths/tests/test_dense.py                                      %{buildroot}/DSR/networkx/algorithms/shortest_paths/tests/test_dense.py                                                                                     
install -m 644 networkx/algorithms/shortest_paths/tests/test_unweighted.py                                 %{buildroot}/DSR/networkx/algorithms/shortest_paths/tests/test_unweighted.py                                                                                
install -m 644 networkx/algorithms/shortest_paths/dense.py                                                 %{buildroot}/DSR/networkx/algorithms/shortest_paths/dense.py                                                                                                
install -m 644 networkx/algorithms/shortest_paths/unweighted.py                                            %{buildroot}/DSR/networkx/algorithms/shortest_paths/unweighted.py                                                                                           
install -m 644 networkx/algorithms/shortest_paths/weighted.py                                              %{buildroot}/DSR/networkx/algorithms/shortest_paths/weighted.py                                                                                             
install -m 644 networkx/algorithms/shortest_paths/astar.py                                                 %{buildroot}/DSR/networkx/algorithms/shortest_paths/astar.py                                                                                                
install -m 644 networkx/algorithms/shortest_paths/__init__.py                                              %{buildroot}/DSR/networkx/algorithms/shortest_paths/__init__.py                                                                                             
install -m 644 networkx/algorithms/shortest_paths/generic.py                                               %{buildroot}/DSR/networkx/algorithms/shortest_paths/generic.py                                                                                              
install -m 644 networkx/algorithms/simple_paths.py                                                         %{buildroot}/DSR/networkx/algorithms/simple_paths.py                                                                                                        
install -m 644 networkx/algorithms/vitality.py                                                             %{buildroot}/DSR/networkx/algorithms/vitality.py                                                                                                            
install -m 644 networkx/algorithms/tests/test_mis.py                                                       %{buildroot}/DSR/networkx/algorithms/tests/test_mis.py                                                                                                      
install -m 644 networkx/algorithms/tests/test_hierarchy.py                                                 %{buildroot}/DSR/networkx/algorithms/tests/test_hierarchy.py                                                                                                
install -m 644 networkx/algorithms/tests/test_swap.py                                                      %{buildroot}/DSR/networkx/algorithms/tests/test_swap.py                                                                                                     
install -m 644 networkx/algorithms/tests/test_simple_paths.py                                              %{buildroot}/DSR/networkx/algorithms/tests/test_simple_paths.py                                                                                             
install -m 644 networkx/algorithms/tests/test_cluster.py                                                   %{buildroot}/DSR/networkx/algorithms/tests/test_cluster.py                                                                                                  
install -m 644 networkx/algorithms/tests/test_boundary.py                                                  %{buildroot}/DSR/networkx/algorithms/tests/test_boundary.py                                                                                                 
install -m 644 networkx/algorithms/tests/test_distance_regular.py                                          %{buildroot}/DSR/networkx/algorithms/tests/test_distance_regular.py                                                                                         
install -m 644 networkx/algorithms/tests/test_euler.py                                                     %{buildroot}/DSR/networkx/algorithms/tests/test_euler.py                                                                                                    
install -m 644 networkx/algorithms/tests/test_graphical.py                                                 %{buildroot}/DSR/networkx/algorithms/tests/test_graphical.py
install -m 644 networkx/algorithms/tests/test_dag.py                                                       %{buildroot}/DSR/networkx/algorithms/tests/test_dag.py
install -m 644 networkx/algorithms/tests/test_smetric.py                                                   %{buildroot}/DSR/networkx/algorithms/tests/test_smetric.py
install -m 644 networkx/algorithms/tests/test_richclub.py                                                  %{buildroot}/DSR/networkx/algorithms/tests/test_richclub.py
install -m 644 networkx/algorithms/tests/test_distance_measures.py                                         %{buildroot}/DSR/networkx/algorithms/tests/test_distance_measures.py
install -m 644 networkx/algorithms/tests/test_matching.py                                                  %{buildroot}/DSR/networkx/algorithms/tests/test_matching.py
install -m 644 networkx/algorithms/tests/test_block.py                                                     %{buildroot}/DSR/networkx/algorithms/tests/test_block.py
install -m 644 networkx/algorithms/tests/test_cycles.py                                                    %{buildroot}/DSR/networkx/algorithms/tests/test_cycles.py
install -m 644 networkx/algorithms/tests/test_clique.py                                                    %{buildroot}/DSR/networkx/algorithms/tests/test_clique.py
install -m 644 networkx/algorithms/tests/test_vitality.py                                                  %{buildroot}/DSR/networkx/algorithms/tests/test_vitality.py
install -m 644 networkx/algorithms/tests/test_mst.py                                                       %{buildroot}/DSR/networkx/algorithms/tests/test_mst.py
install -m 644 networkx/algorithms/tests/test_core.py                                                      %{buildroot}/DSR/networkx/algorithms/tests/test_core.py
install -m 644 networkx/algorithms/approximation/independent_set.py                                        %{buildroot}/DSR/networkx/algorithms/approximation/independent_set.py
install -m 644 networkx/algorithms/approximation/tests/test_independent_set.py                             %{buildroot}/DSR/networkx/algorithms/approximation/tests/test_independent_set.py
install -m 644 networkx/algorithms/approximation/tests/test_vertex_cover.py                                %{buildroot}/DSR/networkx/algorithms/approximation/tests/test_vertex_cover.py
install -m 644 networkx/algorithms/approximation/tests/test_ramsey.py                                      %{buildroot}/DSR/networkx/algorithms/approximation/tests/test_ramsey.py
install -m 644 networkx/algorithms/approximation/tests/test_dominating_set.py                              %{buildroot}/DSR/networkx/algorithms/approximation/tests/test_dominating_set.py
install -m 644 networkx/algorithms/approximation/tests/test_matching.py                                    %{buildroot}/DSR/networkx/algorithms/approximation/tests/test_matching.py
install -m 644 networkx/algorithms/approximation/tests/test_clique.py                                      %{buildroot}/DSR/networkx/algorithms/approximation/tests/test_clique.py
install -m 644 networkx/algorithms/approximation/dominating_set.py                                         %{buildroot}/DSR/networkx/algorithms/approximation/dominating_set.py
install -m 644 networkx/algorithms/approximation/matching.py                                               %{buildroot}/DSR/networkx/algorithms/approximation/matching.py
install -m 644 networkx/algorithms/approximation/ramsey.py                                                 %{buildroot}/DSR/networkx/algorithms/approximation/ramsey.py
install -m 644 networkx/algorithms/approximation/vertex_cover.py                                           %{buildroot}/DSR/networkx/algorithms/approximation/vertex_cover.py
install -m 644 networkx/algorithms/approximation/clique.py                                                 %{buildroot}/DSR/networkx/algorithms/approximation/clique.py
install -m 644 networkx/algorithms/approximation/__init__.py                                               %{buildroot}/DSR/networkx/algorithms/approximation/__init__.py
install -m 644 networkx/algorithms/chordal/tests/test_chordal.py                                           %{buildroot}/DSR/networkx/algorithms/chordal/tests/test_chordal.py
install -m 644 networkx/algorithms/chordal/chordal_alg.py                                                  %{buildroot}/DSR/networkx/algorithms/chordal/chordal_alg.py
install -m 644 networkx/algorithms/chordal/__init__.py                                                     %{buildroot}/DSR/networkx/algorithms/chordal/__init__.py
install -m 644 networkx/algorithms/hierarchy.py                                                            %{buildroot}/DSR/networkx/algorithms/hierarchy.py
install -m 644 networkx/algorithms/block.py                                                                %{buildroot}/DSR/networkx/algorithms/block.py
install -m 644 networkx/algorithms/core.py                                                                 %{buildroot}/DSR/networkx/algorithms/core.py
install -m 644 networkx/algorithms/distance_regular.py                                                     %{buildroot}/DSR/networkx/algorithms/distance_regular.py
install -m 644 networkx/algorithms/components/tests/test_biconnected.py                                    %{buildroot}/DSR/networkx/algorithms/components/tests/test_biconnected.py
install -m 644 networkx/algorithms/components/tests/test_weakly_connected.py                               %{buildroot}/DSR/networkx/algorithms/components/tests/test_weakly_connected.py
install -m 644 networkx/algorithms/components/tests/test_attracting.py                                     %{buildroot}/DSR/networkx/algorithms/components/tests/test_attracting.py
install -m 644 networkx/algorithms/components/tests/test_connected.py                                      %{buildroot}/DSR/networkx/algorithms/components/tests/test_connected.py
install -m 644 networkx/algorithms/components/tests/test_strongly_connected.py                             %{buildroot}/DSR/networkx/algorithms/components/tests/test_strongly_connected.py
install -m 644 networkx/algorithms/components/strongly_connected.py                                        %{buildroot}/DSR/networkx/algorithms/components/strongly_connected.py
install -m 644 networkx/algorithms/components/weakly_connected.py                                          %{buildroot}/DSR/networkx/algorithms/components/weakly_connected.py
install -m 644 networkx/algorithms/components/attracting.py                                                %{buildroot}/DSR/networkx/algorithms/components/attracting.py
install -m 644 networkx/algorithms/components/biconnected.py                                               %{buildroot}/DSR/networkx/algorithms/components/biconnected.py
install -m 644 networkx/algorithms/components/__init__.py                                                  %{buildroot}/DSR/networkx/algorithms/components/__init__.py
install -m 644 networkx/algorithms/components/connected.py                                                 %{buildroot}/DSR/networkx/algorithms/components/connected.py
install -m 644 networkx/algorithms/centrality/closeness.py                                                 %{buildroot}/DSR/networkx/algorithms/centrality/closeness.py
install -m 644 networkx/algorithms/centrality/tests/test_load_centrality.py                                %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_load_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_betweenness_centrality_subset.py                  %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_betweenness_centrality_subset.py
install -m 644 networkx/algorithms/centrality/tests/test_communicability.py                                %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_communicability.py
install -m 644 networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality_subset.py     %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality_subset.py
install -m 644 networkx/algorithms/centrality/tests/test_katz_centrality.py                                %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_katz_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_current_flow_closeness.py                         %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_current_flow_closeness.py
install -m 644 networkx/algorithms/centrality/tests/test_eigenvector_centrality.py                         %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_eigenvector_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality.py            %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_closeness_centrality.py                           %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_closeness_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_betweenness_centrality.py                         %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_betweenness_centrality.py
install -m 644 networkx/algorithms/centrality/tests/test_degree_centrality.py                              %{buildroot}/DSR/networkx/algorithms/centrality/tests/test_degree_centrality.py
install -m 644 networkx/algorithms/centrality/katz.py                                                      %{buildroot}/DSR/networkx/algorithms/centrality/katz.py
install -m 644 networkx/algorithms/centrality/betweenness.py                                               %{buildroot}/DSR/networkx/algorithms/centrality/betweenness.py
install -m 644 networkx/algorithms/centrality/flow_matrix.py                                               %{buildroot}/DSR/networkx/algorithms/centrality/flow_matrix.py
install -m 644 networkx/algorithms/centrality/load.py                                                      %{buildroot}/DSR/networkx/algorithms/centrality/load.py
install -m 644 networkx/algorithms/centrality/current_flow_closeness.py                                    %{buildroot}/DSR/networkx/algorithms/centrality/current_flow_closeness.py
install -m 644 networkx/algorithms/centrality/communicability_alg.py                                       %{buildroot}/DSR/networkx/algorithms/centrality/communicability_alg.py
install -m 644 networkx/algorithms/centrality/eigenvector.py                                               %{buildroot}/DSR/networkx/algorithms/centrality/eigenvector.py
install -m 644 networkx/algorithms/centrality/degree_alg.py                                                %{buildroot}/DSR/networkx/algorithms/centrality/degree_alg.py
install -m 644 networkx/algorithms/centrality/current_flow_betweenness_subset.py                           %{buildroot}/DSR/networkx/algorithms/centrality/current_flow_betweenness_subset.py
install -m 644 networkx/algorithms/centrality/current_flow_betweenness.py                                  %{buildroot}/DSR/networkx/algorithms/centrality/current_flow_betweenness.py
install -m 644 networkx/algorithms/centrality/betweenness_subset.py                                        %{buildroot}/DSR/networkx/algorithms/centrality/betweenness_subset.py
install -m 644 networkx/algorithms/centrality/__init__.py                                                  %{buildroot}/DSR/networkx/algorithms/centrality/__init__.py
install -m 644 networkx/algorithms/matching.py                                                             %{buildroot}/DSR/networkx/algorithms/matching.py
install -m 644 networkx/algorithms/dag.py                                                                  %{buildroot}/DSR/networkx/algorithms/dag.py
install -m 644 networkx/algorithms/boundary.py                                                             %{buildroot}/DSR/networkx/algorithms/boundary.py
install -m 644 networkx/algorithms/euler.py                                                                %{buildroot}/DSR/networkx/algorithms/euler.py
install -m 644 networkx/algorithms/flow/tests/test_mincost.py                                              %{buildroot}/DSR/networkx/algorithms/flow/tests/test_mincost.py
install -m 644 networkx/algorithms/flow/tests/test_maxflow.py                                              %{buildroot}/DSR/networkx/algorithms/flow/tests/test_maxflow.py
install -m 644 networkx/algorithms/flow/tests/test_maxflow_large_graph.py                                  %{buildroot}/DSR/networkx/algorithms/flow/tests/test_maxflow_large_graph.py
install -m 644 networkx/algorithms/flow/mincost.py                                                         %{buildroot}/DSR/networkx/algorithms/flow/mincost.py
install -m 644 networkx/algorithms/flow/maxflow.py                                                         %{buildroot}/DSR/networkx/algorithms/flow/maxflow.py
install -m 644 networkx/algorithms/flow/__init__.py                                                        %{buildroot}/DSR/networkx/algorithms/flow/__init__.py
install -m 644 networkx/algorithms/distance_measures.py                                                    %{buildroot}/DSR/networkx/algorithms/distance_measures.py
install -m 644 networkx/algorithms/isolate.py                                                              %{buildroot}/DSR/networkx/algorithms/isolate.py
install -m 644 networkx/algorithms/community/tests/test_kclique.py                                         %{buildroot}/DSR/networkx/algorithms/community/tests/test_kclique.py
install -m 644 networkx/algorithms/community/kclique.py                                                    %{buildroot}/DSR/networkx/algorithms/community/kclique.py
install -m 644 networkx/algorithms/community/__init__.py                                                   %{buildroot}/DSR/networkx/algorithms/community/__init__.py
install -m 644 networkx/algorithms/swap.py                                                                 %{buildroot}/DSR/networkx/algorithms/swap.py
install -m 644 networkx/algorithms/cycles.py                                                               %{buildroot}/DSR/networkx/algorithms/cycles.py
install -m 644 networkx/algorithms/connectivity/tests/test_cuts.py                                         %{buildroot}/DSR/networkx/algorithms/connectivity/tests/test_cuts.py
install -m 644 networkx/algorithms/connectivity/tests/test_connectivity.py                                 %{buildroot}/DSR/networkx/algorithms/connectivity/tests/test_connectivity.py
install -m 644 networkx/algorithms/connectivity/__init__.py                                                %{buildroot}/DSR/networkx/algorithms/connectivity/__init__.py
install -m 644 networkx/algorithms/connectivity/connectivity.py                                            %{buildroot}/DSR/networkx/algorithms/connectivity/connectivity.py
install -m 644 networkx/algorithms/connectivity/cuts.py                                                    %{buildroot}/DSR/networkx/algorithms/connectivity/cuts.py
install -m 644 networkx/algorithms/clique.py                                                               %{buildroot}/DSR/networkx/algorithms/clique.py
install -m 644 networkx/algorithms/__init__.py                                                             %{buildroot}/DSR/networkx/algorithms/__init__.py
install -m 644 networkx/algorithms/mst.py                                                                  %{buildroot}/DSR/networkx/algorithms/mst.py
install -m 644 networkx/algorithms/bipartite/centrality.py                                                 %{buildroot}/DSR/networkx/algorithms/bipartite/centrality.py
install -m 644 networkx/algorithms/bipartite/projection.py                                                 %{buildroot}/DSR/networkx/algorithms/bipartite/projection.py
install -m 644 networkx/algorithms/bipartite/tests/test_basic.py                                           %{buildroot}/DSR/networkx/algorithms/bipartite/tests/test_basic.py
install -m 644 networkx/algorithms/bipartite/tests/test_centrality.py                                      %{buildroot}/DSR/networkx/algorithms/bipartite/tests/test_centrality.py
install -m 644 networkx/algorithms/bipartite/tests/test_cluster.py                                         %{buildroot}/DSR/networkx/algorithms/bipartite/tests/test_cluster.py
install -m 644 networkx/algorithms/bipartite/tests/test_spectral_bipartivity.py                            %{buildroot}/DSR/networkx/algorithms/bipartite/tests/test_spectral_bipartivity.py
install -m 644 networkx/algorithms/bipartite/tests/test_project.py                                         %{buildroot}/DSR/networkx/algorithms/bipartite/tests/test_project.py
install -m 644 networkx/algorithms/bipartite/basic.py                                                      %{buildroot}/DSR/networkx/algorithms/bipartite/basic.py
install -m 644 networkx/algorithms/bipartite/spectral.py                                                   %{buildroot}/DSR/networkx/algorithms/bipartite/spectral.py
install -m 644 networkx/algorithms/bipartite/__init__.py                                                   %{buildroot}/DSR/networkx/algorithms/bipartite/__init__.py
install -m 644 networkx/algorithms/bipartite/redundancy.py                                                 %{buildroot}/DSR/networkx/algorithms/bipartite/redundancy.py
install -m 644 networkx/algorithms/bipartite/cluster.py                                                    %{buildroot}/DSR/networkx/algorithms/bipartite/cluster.py
install -m 644 networkx/algorithms/traversal/tests/test_bfs.py                                             %{buildroot}/DSR/networkx/algorithms/traversal/tests/test_bfs.py
install -m 644 networkx/algorithms/traversal/tests/test_dfs.py                                             %{buildroot}/DSR/networkx/algorithms/traversal/tests/test_dfs.py
install -m 644 networkx/algorithms/traversal/depth_first_search.py                                         %{buildroot}/DSR/networkx/algorithms/traversal/depth_first_search.py
install -m 644 networkx/algorithms/traversal/__init__.py                                                   %{buildroot}/DSR/networkx/algorithms/traversal/__init__.py
install -m 644 networkx/algorithms/traversal/breadth_first_search.py                                       %{buildroot}/DSR/networkx/algorithms/traversal/breadth_first_search.py
install -m 644 networkx/algorithms/mis.py                                                                  %{buildroot}/DSR/networkx/algorithms/mis.py
install -m 644 networkx/algorithms/operators/unary.py                                                      %{buildroot}/DSR/networkx/algorithms/operators/unary.py
install -m 644 networkx/algorithms/operators/tests/test_product.py                                         %{buildroot}/DSR/networkx/algorithms/operators/tests/test_product.py
install -m 644 networkx/algorithms/operators/tests/test_binary.py                                          %{buildroot}/DSR/networkx/algorithms/operators/tests/test_binary.py
install -m 644 networkx/algorithms/operators/tests/test_all.py                                             %{buildroot}/DSR/networkx/algorithms/operators/tests/test_all.py
install -m 644 networkx/algorithms/operators/tests/test_unary.py                                           %{buildroot}/DSR/networkx/algorithms/operators/tests/test_unary.py
install -m 644 networkx/algorithms/operators/all.py                                                        %{buildroot}/DSR/networkx/algorithms/operators/all.py
install -m 644 networkx/algorithms/operators/product.py                                                    %{buildroot}/DSR/networkx/algorithms/operators/product.py
install -m 644 networkx/algorithms/operators/binary.py                                                     %{buildroot}/DSR/networkx/algorithms/operators/binary.py
install -m 644 networkx/algorithms/operators/__init__.py                                                   %{buildroot}/DSR/networkx/algorithms/operators/__init__.py
install -m 644 networkx/algorithms/isomorphism/isomorphvf2.py                                              %{buildroot}/DSR/networkx/algorithms/isomorphism/isomorphvf2.py
install -m 644 networkx/algorithms/isomorphism/tests/iso_r01_s80.A99                                       %{buildroot}/DSR/networkx/algorithms/isomorphism/tests/iso_r01_s80.A99
install -m 644 networkx/algorithms/isomorphism/tests/test_isomorphvf2.py                                   %{buildroot}/DSR/networkx/algorithms/isomorphism/tests/test_isomorphvf2.py
install -m 644 networkx/algorithms/isomorphism/tests/iso_r01_s80.B99                                       %{buildroot}/DSR/networkx/algorithms/isomorphism/tests/iso_r01_s80.B99
install -m 644 networkx/algorithms/isomorphism/tests/si2_b06_m200.B99                                      %{buildroot}/DSR/networkx/algorithms/isomorphism/tests/si2_b06_m200.B99
install -m 644 networkx/algorithms/isomorphism/tests/test_vf2userfunc.py                                   %{buildroot}/DSR/networkx/algorithms/isomorphism/tests/test_vf2userfunc.py
install -m 644 networkx/algorithms/isomorphism/tests/si2_b06_m200.A99                                      %{buildroot}/DSR/networkx/algorithms/isomorphism/tests/si2_b06_m200.A99
install -m 644 networkx/algorithms/isomorphism/tests/test_isomorphism.py                                   %{buildroot}/DSR/networkx/algorithms/isomorphism/tests/test_isomorphism.py
install -m 644 networkx/algorithms/isomorphism/matchhelpers.py                                             %{buildroot}/DSR/networkx/algorithms/isomorphism/matchhelpers.py
install -m 644 networkx/algorithms/isomorphism/__init__.py                                                 %{buildroot}/DSR/networkx/algorithms/isomorphism/__init__.py
install -m 644 networkx/algorithms/isomorphism/isomorph.py                                                 %{buildroot}/DSR/networkx/algorithms/isomorphism/isomorph.py
install -m 644 networkx/algorithms/isomorphism/vf2userfunc.py                                              %{buildroot}/DSR/networkx/algorithms/isomorphism/vf2userfunc.py
install -m 644 networkx/algorithms/cluster.py                                                              %{buildroot}/DSR/networkx/algorithms/cluster.py
install -m 644 networkx/tests/test_convert.py                                                              %{buildroot}/DSR/networkx/tests/test_convert.py
install -m 644 networkx/tests/test_convert_numpy.py                                                        %{buildroot}/DSR/networkx/tests/test_convert_numpy.py
install -m 644 networkx/tests/test_exceptions.py                                                           %{buildroot}/DSR/networkx/tests/test_exceptions.py
install -m 644 networkx/tests/test_relabel.py                                                              %{buildroot}/DSR/networkx/tests/test_relabel.py
install -m 644 networkx/tests/test.py                                                                      %{buildroot}/DSR/networkx/tests/test.py
install -m 644 networkx/tests/__init__.py                                                                  %{buildroot}/DSR/networkx/tests/__init__.py
install -m 644 networkx/tests/test_convert_scipy.py                                                        %{buildroot}/DSR/networkx/tests/test_convert_scipy.py
install -m 644 networkx/tests/benchmark.py                                                                 %{buildroot}/DSR/networkx/tests/benchmark.py
install -m 644 networkx/generators/tests/test_random_graphs.py                                             %{buildroot}/DSR/networkx/generators/tests/test_random_graphs.py
install -m 644 networkx/generators/tests/test_small.py                                                     %{buildroot}/DSR/networkx/generators/tests/test_small.py
install -m 644 networkx/generators/tests/test_geometric.py                                                 %{buildroot}/DSR/networkx/generators/tests/test_geometric.py
install -m 644 networkx/generators/tests/test_stochastic.py                                                %{buildroot}/DSR/networkx/generators/tests/test_stochastic.py
install -m 644 networkx/generators/tests/test_classic.py                                                   %{buildroot}/DSR/networkx/generators/tests/test_classic.py
install -m 644 networkx/generators/tests/test_threshold.py                                                 %{buildroot}/DSR/networkx/generators/tests/test_threshold.py
install -m 644 networkx/generators/tests/test_atlas.py                                                     %{buildroot}/DSR/networkx/generators/tests/test_atlas.py
install -m 644 networkx/generators/tests/test_line.py                                                      %{buildroot}/DSR/networkx/generators/tests/test_line.py
install -m 644 networkx/generators/tests/test_random_clustered.py                                          %{buildroot}/DSR/networkx/generators/tests/test_random_clustered.py
install -m 644 networkx/generators/tests/test_hybrid.py                                                    %{buildroot}/DSR/networkx/generators/tests/test_hybrid.py
install -m 644 networkx/generators/tests/test_directed.py                                                  %{buildroot}/DSR/networkx/generators/tests/test_directed.py
install -m 644 networkx/generators/tests/test_degree_seq.py                                                %{buildroot}/DSR/networkx/generators/tests/test_degree_seq.py
install -m 644 networkx/generators/tests/test_bipartite.py                                                 %{buildroot}/DSR/networkx/generators/tests/test_bipartite.py
install -m 644 networkx/generators/tests/test_intersection.py                                              %{buildroot}/DSR/networkx/generators/tests/test_intersection.py
install -m 644 networkx/generators/tests/test_ego.py                                                       %{buildroot}/DSR/networkx/generators/tests/test_ego.py
install -m 644 networkx/generators/classic.py                                                              %{buildroot}/DSR/networkx/generators/classic.py
install -m 644 networkx/generators/social.py                                                               %{buildroot}/DSR/networkx/generators/social.py
install -m 644 networkx/generators/intersection.py                                                         %{buildroot}/DSR/networkx/generators/intersection.py
install -m 644 networkx/generators/threshold.py                                                            %{buildroot}/DSR/networkx/generators/threshold.py
install -m 644 networkx/generators/ego.py                                                                  %{buildroot}/DSR/networkx/generators/ego.py
install -m 644 networkx/generators/atlas.py                                                                %{buildroot}/DSR/networkx/generators/atlas.py
install -m 644 networkx/generators/hybrid.py                                                               %{buildroot}/DSR/networkx/generators/hybrid.py
install -m 644 networkx/generators/directed.py                                                             %{buildroot}/DSR/networkx/generators/directed.py
install -m 644 networkx/generators/random_graphs.py                                                        %{buildroot}/DSR/networkx/generators/random_graphs.py
install -m 644 networkx/generators/bipartite.py                                                            %{buildroot}/DSR/networkx/generators/bipartite.py
install -m 644 networkx/generators/line.py                                                                 %{buildroot}/DSR/networkx/generators/line.py
install -m 644 networkx/generators/__init__.py                                                             %{buildroot}/DSR/networkx/generators/__init__.py
install -m 644 networkx/generators/geometric.py                                                            %{buildroot}/DSR/networkx/generators/geometric.py
install -m 644 networkx/generators/random_clustered.py                                                     %{buildroot}/DSR/networkx/generators/random_clustered.py
install -m 644 networkx/generators/small.py                                                                %{buildroot}/DSR/networkx/generators/small.py
install -m 644 networkx/generators/degree_seq.py                                                           %{buildroot}/DSR/networkx/generators/degree_seq.py
install -m 644 networkx/generators/stochastic.py                                                           %{buildroot}/DSR/networkx/generators/stochastic.py
install -m 644 networkx/exception.py                                                                       %{buildroot}/DSR/networkx/exception.py
install -m 644 networkx/version.py                                                                         %{buildroot}/DSR/networkx/version.py
install -m 644 networkx/release.py                                                                         %{buildroot}/DSR/networkx/release.py
install -m 644 networkx/linalg/tests/test_spectrum.py                                                      %{buildroot}/DSR/networkx/linalg/tests/test_spectrum.py
install -m 644 networkx/linalg/tests/test_laplacian.py                                                     %{buildroot}/DSR/networkx/linalg/tests/test_laplacian.py
install -m 644 networkx/linalg/tests/test_graphmatrix.py                                                   %{buildroot}/DSR/networkx/linalg/tests/test_graphmatrix.py
install -m 644 networkx/linalg/attrmatrix.py                                                               %{buildroot}/DSR/networkx/linalg/attrmatrix.py
install -m 644 networkx/linalg/graphmatrix.py                                                              %{buildroot}/DSR/networkx/linalg/graphmatrix.py
install -m 644 networkx/linalg/__init__.py                                                                 %{buildroot}/DSR/networkx/linalg/__init__.py
install -m 644 networkx/linalg/laplacianmatrix.py                                                          %{buildroot}/DSR/networkx/linalg/laplacianmatrix.py
install -m 644 networkx/linalg/spectrum.py                                                                 %{buildroot}/DSR/networkx/linalg/spectrum.py
install -m 644 networkx/readwrite/leda.py                                                                  %{buildroot}/DSR/networkx/readwrite/leda.py
install -m 644 networkx/readwrite/tests/test_p2g.py                                                        %{buildroot}/DSR/networkx/readwrite/tests/test_p2g.py
install -m 644 networkx/readwrite/tests/test_leda.py                                                       %{buildroot}/DSR/networkx/readwrite/tests/test_leda.py
install -m 644 networkx/readwrite/tests/test_shp.py                                                        %{buildroot}/DSR/networkx/readwrite/tests/test_shp.py
install -m 644 networkx/readwrite/tests/test_gpickle.py                                                    %{buildroot}/DSR/networkx/readwrite/tests/test_gpickle.py
install -m 644 networkx/readwrite/tests/test_yaml.py                                                       %{buildroot}/DSR/networkx/readwrite/tests/test_yaml.py
install -m 644 networkx/readwrite/tests/test_graphml.py                                                    %{buildroot}/DSR/networkx/readwrite/tests/test_graphml.py
install -m 644 networkx/readwrite/tests/test_gexf.py                                                       %{buildroot}/DSR/networkx/readwrite/tests/test_gexf.py
install -m 644 networkx/readwrite/tests/test_sparsegraph6.py                                               %{buildroot}/DSR/networkx/readwrite/tests/test_sparsegraph6.py
install -m 644 networkx/readwrite/tests/test_gml.py                                                        %{buildroot}/DSR/networkx/readwrite/tests/test_gml.py
install -m 644 networkx/readwrite/tests/test_adjlist.py                                                    %{buildroot}/DSR/networkx/readwrite/tests/test_adjlist.py
install -m 644 networkx/readwrite/tests/test_edgelist.py                                                   %{buildroot}/DSR/networkx/readwrite/tests/test_edgelist.py
install -m 644 networkx/readwrite/tests/test_pajek.py                                                      %{buildroot}/DSR/networkx/readwrite/tests/test_pajek.py
install -m 644 networkx/readwrite/multiline_adjlist.py                                                     %{buildroot}/DSR/networkx/readwrite/multiline_adjlist.py
install -m 644 networkx/readwrite/json_graph/tests/test_tree.py                                            %{buildroot}/DSR/networkx/readwrite/json_graph/tests/test_tree.py
install -m 644 networkx/readwrite/json_graph/tests/test_adjacency.py                                       %{buildroot}/DSR/networkx/readwrite/json_graph/tests/test_adjacency.py
install -m 644 networkx/readwrite/json_graph/tests/test_serialize.py                                       %{buildroot}/DSR/networkx/readwrite/json_graph/tests/test_serialize.py
install -m 644 networkx/readwrite/json_graph/tests/test_node_link.py                                       %{buildroot}/DSR/networkx/readwrite/json_graph/tests/test_node_link.py
install -m 644 networkx/readwrite/json_graph/tree.py                                                       %{buildroot}/DSR/networkx/readwrite/json_graph/tree.py
install -m 644 networkx/readwrite/json_graph/serialize.py                                                  %{buildroot}/DSR/networkx/readwrite/json_graph/serialize.py
install -m 644 networkx/readwrite/json_graph/node_link.py                                                  %{buildroot}/DSR/networkx/readwrite/json_graph/node_link.py
install -m 644 networkx/readwrite/json_graph/__init__.py                                                   %{buildroot}/DSR/networkx/readwrite/json_graph/__init__.py
install -m 644 networkx/readwrite/json_graph/adjacency.py                                                  %{buildroot}/DSR/networkx/readwrite/json_graph/adjacency.py
install -m 644 networkx/readwrite/edgelist.py                                                              %{buildroot}/DSR/networkx/readwrite/edgelist.py
install -m 644 networkx/readwrite/gml.py                                                                   %{buildroot}/DSR/networkx/readwrite/gml.py
install -m 644 networkx/readwrite/pajek.py                                                                 %{buildroot}/DSR/networkx/readwrite/pajek.py
install -m 644 networkx/readwrite/graphml.py                                                               %{buildroot}/DSR/networkx/readwrite/graphml.py
install -m 644 networkx/readwrite/gexf.py                                                                  %{buildroot}/DSR/networkx/readwrite/gexf.py
install -m 644 networkx/readwrite/adjlist.py                                                               %{buildroot}/DSR/networkx/readwrite/adjlist.py
install -m 644 networkx/readwrite/gpickle.py                                                               %{buildroot}/DSR/networkx/readwrite/gpickle.py
install -m 644 networkx/readwrite/__init__.py                                                              %{buildroot}/DSR/networkx/readwrite/__init__.py
install -m 644 networkx/readwrite/sparsegraph6.py                                                          %{buildroot}/DSR/networkx/readwrite/sparsegraph6.py
install -m 644 networkx/readwrite/nx_yaml.py                                                               %{buildroot}/DSR/networkx/readwrite/nx_yaml.py
install -m 644 networkx/readwrite/p2g.py                                                                   %{buildroot}/DSR/networkx/readwrite/p2g.py
install -m 644 networkx/readwrite/nx_shp.py                                                                %{buildroot}/DSR/networkx/readwrite/nx_shp.py
install -m 644 networkx/relabel.py                                                                         %{buildroot}/DSR/networkx/relabel.py
install -m 644 networkx/classes/tests/test_multigraph.py                                                   %{buildroot}/DSR/networkx/classes/tests/test_multigraph.py
install -m 644 networkx/classes/tests/test_graph.py                                                        %{buildroot}/DSR/networkx/classes/tests/test_graph.py
install -m 644 networkx/classes/tests/test_multidigraph.py                                                 %{buildroot}/DSR/networkx/classes/tests/test_multidigraph.py
install -m 644 networkx/classes/tests/test_digraph.py                                                      %{buildroot}/DSR/networkx/classes/tests/test_digraph.py
install -m 644 networkx/classes/tests/test_graph_historical.py                                             %{buildroot}/DSR/networkx/classes/tests/test_graph_historical.py
install -m 644 networkx/classes/tests/historical_tests.py                                                  %{buildroot}/DSR/networkx/classes/tests/historical_tests.py
install -m 644 networkx/classes/tests/test_digraph_historical.py                                           %{buildroot}/DSR/networkx/classes/tests/test_digraph_historical.py
install -m 644 networkx/classes/tests/test_function.py                                                     %{buildroot}/DSR/networkx/classes/tests/test_function.py
install -m 644 networkx/classes/graph.py                                                                   %{buildroot}/DSR/networkx/classes/graph.py
install -m 644 networkx/classes/multigraph.py                                                              %{buildroot}/DSR/networkx/classes/multigraph.py
install -m 644 networkx/classes/function.py                                                                %{buildroot}/DSR/networkx/classes/function.py
install -m 644 networkx/classes/digraph.py                                                                 %{buildroot}/DSR/networkx/classes/digraph.py
install -m 644 networkx/classes/__init__.py                                                                %{buildroot}/DSR/networkx/classes/__init__.py
install -m 644 networkx/classes/multidigraph.py                                                            %{buildroot}/DSR/networkx/classes/multidigraph.py
install -m 644 networkx/utils/tests/test_misc.py                                                           %{buildroot}/DSR/networkx/utils/tests/test_misc.py
install -m 644 networkx/utils/tests/test_decorators.py                                                     %{buildroot}/DSR/networkx/utils/tests/test_decorators.py
install -m 644 networkx/utils/tests/test_random_sequence.py                                                %{buildroot}/DSR/networkx/utils/tests/test_random_sequence.py
install -m 644 networkx/utils/tests/test_rcm.py                                                            %{buildroot}/DSR/networkx/utils/tests/test_rcm.py
install -m 644 networkx/utils/tests/test.txt                                                               %{buildroot}/DSR/networkx/utils/tests/test.txt
install -m 644 networkx/utils/union_find.py                                                                %{buildroot}/DSR/networkx/utils/union_find.py
install -m 644 networkx/utils/decorators.py                                                                %{buildroot}/DSR/networkx/utils/decorators.py
install -m 644 networkx/utils/rcm.py                                                                       %{buildroot}/DSR/networkx/utils/rcm.py
install -m 644 networkx/utils/misc.py                                                                      %{buildroot}/DSR/networkx/utils/misc.py
install -m 644 networkx/utils/random_sequence.py                                                           %{buildroot}/DSR/networkx/utils/random_sequence.py
install -m 644 networkx/utils/__init__.py                                                                  %{buildroot}/DSR/networkx/utils/__init__.py
install -m 644 networkx/drawing/tests/test_agraph.py                                                       %{buildroot}/DSR/networkx/drawing/tests/test_agraph.py
install -m 644 networkx/drawing/tests/test_pydot.py                                                        %{buildroot}/DSR/networkx/drawing/tests/test_pydot.py
install -m 644 networkx/drawing/tests/test_layout.py                                                       %{buildroot}/DSR/networkx/drawing/tests/test_layout.py
install -m 644 networkx/drawing/tests/test_pylab.py                                                        %{buildroot}/DSR/networkx/drawing/tests/test_pylab.py
install -m 644 networkx/drawing/nx_agraph.py                                                               %{buildroot}/DSR/networkx/drawing/nx_agraph.py
install -m 644 networkx/drawing/nx_pydot.py                                                                %{buildroot}/DSR/networkx/drawing/nx_pydot.py
install -m 644 networkx/drawing/__init__.py                                                                %{buildroot}/DSR/networkx/drawing/__init__.py
install -m 644 networkx/drawing/layout.py                                                                  %{buildroot}/DSR/networkx/drawing/layout.py
install -m 644 networkx/drawing/nx_pylab.py                                                                %{buildroot}/DSR/networkx/drawing/nx_pylab.py
install -m 644 networkx/__init__.py                                                                        %{buildroot}/DSR/networkx/__init__.py
install -m 644 networkx/convert.py                                                                         %{buildroot}/DSR/networkx/convert.py
install -m 644 networkx/testing/tests/test_utils.py                                                        %{buildroot}/DSR/networkx/testing/tests/test_utils.py
install -m 644 networkx/testing/utils.py                                                                   %{buildroot}/DSR/networkx/testing/utils.py
install -m 644 networkx/testing/__init__.py                                                                %{buildroot}/DSR/networkx/testing/__init__.py

#dos2unix %{buildroot}/*

%files
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
/DSR/example/
/DSR/example/p21c.hkl
/DSR/example/p21c.res
/DSR/example/p21c_step0.res
/DSR/example/p21c_step1.res
/DSR/example/p21c_step2.res
/DSR/example/p21c_step3.res
/DSR/example/p21c-step2.ins

/DSR/networkx/external/decorator/decorator3/_decorator3.py
/DSR/networkx/external/decorator/decorator3/__init__.py
/DSR/networkx/external/decorator/__init__.py
/DSR/networkx/external/decorator/decorator2/_decorator2.py
/DSR/networkx/external/decorator/decorator2/__init__.py
/DSR/networkx/external/__init__.py
/DSR/networkx/algorithms/richclub.py
/DSR/networkx/algorithms/assortativity/tests/test_pairs.py
/DSR/networkx/algorithms/assortativity/tests/test_correlation.py
/DSR/networkx/algorithms/assortativity/tests/test_neighbor_degree.py
/DSR/networkx/algorithms/assortativity/tests/test_mixing.py
/DSR/networkx/algorithms/assortativity/tests/base_test.py
/DSR/networkx/algorithms/assortativity/tests/test_connectivity.py
/DSR/networkx/algorithms/assortativity/pairs.py
/DSR/networkx/algorithms/assortativity/correlation.py
/DSR/networkx/algorithms/assortativity/mixing.py
/DSR/networkx/algorithms/assortativity/__init__.py
/DSR/networkx/algorithms/assortativity/neighbor_degree.py                                                                                            
/DSR/networkx/algorithms/assortativity/connectivity.py                                                                                               
/DSR/networkx/algorithms/link_analysis/tests/test_hits.py                                                                                            
/DSR/networkx/algorithms/link_analysis/tests/test_pagerank.py                                                                                        
/DSR/networkx/algorithms/link_analysis/pagerank_alg.py                                                                                               
/DSR/networkx/algorithms/link_analysis/__init__.py                                                                                                   
/DSR/networkx/algorithms/link_analysis/hits_alg.py                                                                                                   
/DSR/networkx/algorithms/smetric.py                                                                                                                  
/DSR/networkx/algorithms/graphical.py                                                                                                                
/DSR/networkx/algorithms/shortest_paths/tests/test_generic.py                                                                                        
/DSR/networkx/algorithms/shortest_paths/tests/test_weighted.py                                                                                       
/DSR/networkx/algorithms/shortest_paths/tests/test_astar.py                                                                                          
/DSR/networkx/algorithms/shortest_paths/tests/test_dense_numpy.py                                                                                    
/DSR/networkx/algorithms/shortest_paths/tests/test_dense.py                                                                                          
/DSR/networkx/algorithms/shortest_paths/tests/test_unweighted.py                                                                                     
/DSR/networkx/algorithms/shortest_paths/dense.py                                                                                                     
/DSR/networkx/algorithms/shortest_paths/unweighted.py                                                                                                
/DSR/networkx/algorithms/shortest_paths/weighted.py                                                                                                  
/DSR/networkx/algorithms/shortest_paths/astar.py                                                                                                     
/DSR/networkx/algorithms/shortest_paths/__init__.py                                                                                                  
/DSR/networkx/algorithms/shortest_paths/generic.py                                                                                                   
/DSR/networkx/algorithms/simple_paths.py                                                                                                             
/DSR/networkx/algorithms/vitality.py                                                                                                                 
/DSR/networkx/algorithms/tests/test_mis.py                                                                                                           
/DSR/networkx/algorithms/tests/test_hierarchy.py                                                                                                     
/DSR/networkx/algorithms/tests/test_swap.py                                                                                                          
/DSR/networkx/algorithms/tests/test_simple_paths.py                                                                                                  
/DSR/networkx/algorithms/tests/test_cluster.py                                                                                                       
/DSR/networkx/algorithms/tests/test_boundary.py                                                                                                      
/DSR/networkx/algorithms/tests/test_distance_regular.py                                                                                              
/DSR/networkx/algorithms/tests/test_euler.py                                                                                                         
/DSR/networkx/algorithms/tests/test_graphical.py
/DSR/networkx/algorithms/tests/test_dag.py
/DSR/networkx/algorithms/tests/test_smetric.py
/DSR/networkx/algorithms/tests/test_richclub.py
/DSR/networkx/algorithms/tests/test_distance_measures.py
/DSR/networkx/algorithms/tests/test_matching.py
/DSR/networkx/algorithms/tests/test_block.py
/DSR/networkx/algorithms/tests/test_cycles.py
/DSR/networkx/algorithms/tests/test_clique.py
/DSR/networkx/algorithms/tests/test_vitality.py
/DSR/networkx/algorithms/tests/test_mst.py
/DSR/networkx/algorithms/tests/test_core.py
/DSR/networkx/algorithms/approximation/independent_set.py
/DSR/networkx/algorithms/approximation/tests/test_independent_set.py
/DSR/networkx/algorithms/approximation/tests/test_vertex_cover.py
/DSR/networkx/algorithms/approximation/tests/test_ramsey.py
/DSR/networkx/algorithms/approximation/tests/test_dominating_set.py
/DSR/networkx/algorithms/approximation/tests/test_matching.py
/DSR/networkx/algorithms/approximation/tests/test_clique.py
/DSR/networkx/algorithms/approximation/dominating_set.py
/DSR/networkx/algorithms/approximation/matching.py
/DSR/networkx/algorithms/approximation/ramsey.py
/DSR/networkx/algorithms/approximation/vertex_cover.py
/DSR/networkx/algorithms/approximation/clique.py
/DSR/networkx/algorithms/approximation/__init__.py
/DSR/networkx/algorithms/chordal/tests/test_chordal.py
/DSR/networkx/algorithms/chordal/chordal_alg.py
/DSR/networkx/algorithms/chordal/__init__.py
/DSR/networkx/algorithms/hierarchy.py
/DSR/networkx/algorithms/block.py
/DSR/networkx/algorithms/core.py
/DSR/networkx/algorithms/distance_regular.py
/DSR/networkx/algorithms/components/tests/test_biconnected.py
/DSR/networkx/algorithms/components/tests/test_weakly_connected.py
/DSR/networkx/algorithms/components/tests/test_attracting.py
/DSR/networkx/algorithms/components/tests/test_connected.py
/DSR/networkx/algorithms/components/tests/test_strongly_connected.py
/DSR/networkx/algorithms/components/strongly_connected.py
/DSR/networkx/algorithms/components/weakly_connected.py
/DSR/networkx/algorithms/components/attracting.py
/DSR/networkx/algorithms/components/biconnected.py
/DSR/networkx/algorithms/components/__init__.py
/DSR/networkx/algorithms/components/connected.py
/DSR/networkx/algorithms/centrality/closeness.py
/DSR/networkx/algorithms/centrality/tests/test_load_centrality.py
/DSR/networkx/algorithms/centrality/tests/test_betweenness_centrality_subset.py
/DSR/networkx/algorithms/centrality/tests/test_communicability.py
/DSR/networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality_subset.py
/DSR/networkx/algorithms/centrality/tests/test_katz_centrality.py
/DSR/networkx/algorithms/centrality/tests/test_current_flow_closeness.py
/DSR/networkx/algorithms/centrality/tests/test_eigenvector_centrality.py
/DSR/networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality.py
/DSR/networkx/algorithms/centrality/tests/test_closeness_centrality.py
/DSR/networkx/algorithms/centrality/tests/test_betweenness_centrality.py
/DSR/networkx/algorithms/centrality/tests/test_degree_centrality.py
/DSR/networkx/algorithms/centrality/katz.py
/DSR/networkx/algorithms/centrality/betweenness.py
/DSR/networkx/algorithms/centrality/flow_matrix.py
/DSR/networkx/algorithms/centrality/load.py
/DSR/networkx/algorithms/centrality/current_flow_closeness.py
/DSR/networkx/algorithms/centrality/communicability_alg.py
/DSR/networkx/algorithms/centrality/eigenvector.py
/DSR/networkx/algorithms/centrality/degree_alg.py
/DSR/networkx/algorithms/centrality/current_flow_betweenness_subset.py
/DSR/networkx/algorithms/centrality/current_flow_betweenness.py
/DSR/networkx/algorithms/centrality/betweenness_subset.py
/DSR/networkx/algorithms/centrality/__init__.py
/DSR/networkx/algorithms/matching.py
/DSR/networkx/algorithms/dag.py
/DSR/networkx/algorithms/boundary.py
/DSR/networkx/algorithms/euler.py
/DSR/networkx/algorithms/flow/tests/test_mincost.py
/DSR/networkx/algorithms/flow/tests/test_maxflow.py
/DSR/networkx/algorithms/flow/tests/test_maxflow_large_graph.py
/DSR/networkx/algorithms/flow/mincost.py
/DSR/networkx/algorithms/flow/maxflow.py
/DSR/networkx/algorithms/flow/__init__.py
/DSR/networkx/algorithms/distance_measures.py
/DSR/networkx/algorithms/isolate.py
/DSR/networkx/algorithms/community/tests/test_kclique.py
/DSR/networkx/algorithms/community/kclique.py
/DSR/networkx/algorithms/community/__init__.py
/DSR/networkx/algorithms/swap.py
/DSR/networkx/algorithms/cycles.py
/DSR/networkx/algorithms/connectivity/tests/test_cuts.py
/DSR/networkx/algorithms/connectivity/tests/test_connectivity.py
/DSR/networkx/algorithms/connectivity/__init__.py
/DSR/networkx/algorithms/connectivity/connectivity.py
/DSR/networkx/algorithms/connectivity/cuts.py
/DSR/networkx/algorithms/clique.py
/DSR/networkx/algorithms/__init__.py
/DSR/networkx/algorithms/mst.py
/DSR/networkx/algorithms/bipartite/centrality.py
/DSR/networkx/algorithms/bipartite/projection.py
/DSR/networkx/algorithms/bipartite/tests/test_basic.py
/DSR/networkx/algorithms/bipartite/tests/test_centrality.py
/DSR/networkx/algorithms/bipartite/tests/test_cluster.py
/DSR/networkx/algorithms/bipartite/tests/test_spectral_bipartivity.py
/DSR/networkx/algorithms/bipartite/tests/test_project.py
/DSR/networkx/algorithms/bipartite/basic.py
/DSR/networkx/algorithms/bipartite/spectral.py
/DSR/networkx/algorithms/bipartite/__init__.py
/DSR/networkx/algorithms/bipartite/redundancy.py
/DSR/networkx/algorithms/bipartite/cluster.py
/DSR/networkx/algorithms/traversal/tests/test_bfs.py
/DSR/networkx/algorithms/traversal/tests/test_dfs.py
/DSR/networkx/algorithms/traversal/depth_first_search.py
/DSR/networkx/algorithms/traversal/__init__.py
/DSR/networkx/algorithms/traversal/breadth_first_search.py
/DSR/networkx/algorithms/mis.py
/DSR/networkx/algorithms/operators/unary.py
/DSR/networkx/algorithms/operators/tests/test_product.py
/DSR/networkx/algorithms/operators/tests/test_binary.py
/DSR/networkx/algorithms/operators/tests/test_all.py
/DSR/networkx/algorithms/operators/tests/test_unary.py
/DSR/networkx/algorithms/operators/all.py
/DSR/networkx/algorithms/operators/product.py
/DSR/networkx/algorithms/operators/binary.py
/DSR/networkx/algorithms/operators/__init__.py
/DSR/networkx/algorithms/isomorphism/isomorphvf2.py
/DSR/networkx/algorithms/isomorphism/tests/iso_r01_s80.A99
/DSR/networkx/algorithms/isomorphism/tests/test_isomorphvf2.py
/DSR/networkx/algorithms/isomorphism/tests/iso_r01_s80.B99
/DSR/networkx/algorithms/isomorphism/tests/si2_b06_m200.B99
/DSR/networkx/algorithms/isomorphism/tests/test_vf2userfunc.py
/DSR/networkx/algorithms/isomorphism/tests/si2_b06_m200.A99
/DSR/networkx/algorithms/isomorphism/tests/test_isomorphism.py
/DSR/networkx/algorithms/isomorphism/matchhelpers.py
/DSR/networkx/algorithms/isomorphism/__init__.py
/DSR/networkx/algorithms/isomorphism/isomorph.py
/DSR/networkx/algorithms/isomorphism/vf2userfunc.py
/DSR/networkx/algorithms/cluster.py
/DSR/networkx/tests/test_convert.py
/DSR/networkx/tests/test_convert_numpy.py
/DSR/networkx/tests/test_exceptions.py
/DSR/networkx/tests/test_relabel.py
/DSR/networkx/tests/test.py
/DSR/networkx/tests/__init__.py
/DSR/networkx/tests/test_convert_scipy.py
/DSR/networkx/tests/benchmark.py
/DSR/networkx/generators/tests/test_random_graphs.py
/DSR/networkx/generators/tests/test_small.py
/DSR/networkx/generators/tests/test_geometric.py
/DSR/networkx/generators/tests/test_stochastic.py
/DSR/networkx/generators/tests/test_classic.py
/DSR/networkx/generators/tests/test_threshold.py
/DSR/networkx/generators/tests/test_atlas.py
/DSR/networkx/generators/tests/test_line.py
/DSR/networkx/generators/tests/test_random_clustered.py
/DSR/networkx/generators/tests/test_hybrid.py
/DSR/networkx/generators/tests/test_directed.py
/DSR/networkx/generators/tests/test_degree_seq.py
/DSR/networkx/generators/tests/test_bipartite.py
/DSR/networkx/generators/tests/test_intersection.py
/DSR/networkx/generators/tests/test_ego.py
/DSR/networkx/generators/classic.py
/DSR/networkx/generators/social.py
/DSR/networkx/generators/intersection.py
/DSR/networkx/generators/threshold.py
/DSR/networkx/generators/ego.py
/DSR/networkx/generators/atlas.py
/DSR/networkx/generators/hybrid.py
/DSR/networkx/generators/directed.py
/DSR/networkx/generators/random_graphs.py
/DSR/networkx/generators/bipartite.py
/DSR/networkx/generators/line.py
/DSR/networkx/generators/__init__.py
/DSR/networkx/generators/geometric.py
/DSR/networkx/generators/random_clustered.py
/DSR/networkx/generators/small.py
/DSR/networkx/generators/degree_seq.py
/DSR/networkx/generators/stochastic.py
/DSR/networkx/exception.py
/DSR/networkx/version.py
/DSR/networkx/release.py
/DSR/networkx/linalg/tests/test_spectrum.py
/DSR/networkx/linalg/tests/test_laplacian.py
/DSR/networkx/linalg/tests/test_graphmatrix.py
/DSR/networkx/linalg/attrmatrix.py
/DSR/networkx/linalg/graphmatrix.py
/DSR/networkx/linalg/__init__.py
/DSR/networkx/linalg/laplacianmatrix.py
/DSR/networkx/linalg/spectrum.py
/DSR/networkx/readwrite/leda.py
/DSR/networkx/readwrite/tests/test_p2g.py
/DSR/networkx/readwrite/tests/test_leda.py
/DSR/networkx/readwrite/tests/test_shp.py
/DSR/networkx/readwrite/tests/test_gpickle.py
/DSR/networkx/readwrite/tests/test_yaml.py
/DSR/networkx/readwrite/tests/test_graphml.py
/DSR/networkx/readwrite/tests/test_gexf.py
/DSR/networkx/readwrite/tests/test_sparsegraph6.py
/DSR/networkx/readwrite/tests/test_gml.py
/DSR/networkx/readwrite/tests/test_adjlist.py
/DSR/networkx/readwrite/tests/test_edgelist.py
/DSR/networkx/readwrite/tests/test_pajek.py
/DSR/networkx/readwrite/multiline_adjlist.py
/DSR/networkx/readwrite/json_graph/tests/test_tree.py
/DSR/networkx/readwrite/json_graph/tests/test_adjacency.py
/DSR/networkx/readwrite/json_graph/tests/test_serialize.py
/DSR/networkx/readwrite/json_graph/tests/test_node_link.py
/DSR/networkx/readwrite/json_graph/tree.py
/DSR/networkx/readwrite/json_graph/serialize.py
/DSR/networkx/readwrite/json_graph/node_link.py
/DSR/networkx/readwrite/json_graph/__init__.py
/DSR/networkx/readwrite/json_graph/adjacency.py
/DSR/networkx/readwrite/edgelist.py
/DSR/networkx/readwrite/gml.py
/DSR/networkx/readwrite/pajek.py
/DSR/networkx/readwrite/graphml.py
/DSR/networkx/readwrite/gexf.py
/DSR/networkx/readwrite/adjlist.py
/DSR/networkx/readwrite/gpickle.py
/DSR/networkx/readwrite/__init__.py
/DSR/networkx/readwrite/sparsegraph6.py
/DSR/networkx/readwrite/nx_yaml.py
/DSR/networkx/readwrite/p2g.py
/DSR/networkx/readwrite/nx_shp.py
/DSR/networkx/relabel.py
/DSR/networkx/classes/tests/test_multigraph.py
/DSR/networkx/classes/tests/test_graph.py
/DSR/networkx/classes/tests/test_multidigraph.py
/DSR/networkx/classes/tests/test_digraph.py
/DSR/networkx/classes/tests/test_graph_historical.py
/DSR/networkx/classes/tests/historical_tests.py
/DSR/networkx/classes/tests/test_digraph_historical.py
/DSR/networkx/classes/tests/test_function.py
/DSR/networkx/classes/graph.py
/DSR/networkx/classes/multigraph.py
/DSR/networkx/classes/function.py
/DSR/networkx/classes/digraph.py
/DSR/networkx/classes/__init__.py
/DSR/networkx/classes/multidigraph.py
/DSR/networkx/utils/tests/test_misc.py
/DSR/networkx/utils/tests/test_decorators.py
/DSR/networkx/utils/tests/test_random_sequence.py
/DSR/networkx/utils/tests/test_rcm.py
/DSR/networkx/utils/tests/test.txt
/DSR/networkx/utils/union_find.py
/DSR/networkx/utils/decorators.py
/DSR/networkx/utils/rcm.py
/DSR/networkx/utils/misc.py
/DSR/networkx/utils/random_sequence.py
/DSR/networkx/utils/__init__.py
/DSR/networkx/drawing/tests/test_agraph.py
/DSR/networkx/drawing/tests/test_pydot.py
/DSR/networkx/drawing/tests/test_layout.py
/DSR/networkx/drawing/tests/test_pylab.py
/DSR/networkx/drawing/nx_agraph.py
/DSR/networkx/drawing/nx_pydot.py
/DSR/networkx/drawing/__init__.py
/DSR/networkx/drawing/layout.py
/DSR/networkx/drawing/nx_pylab.py
/DSR/networkx/__init__.py
/DSR/networkx/convert.py
/DSR/networkx/testing/tests/test_utils.py
/DSR/networkx/testing/utils.py
/DSR/networkx/testing/__init__.py