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

./networkx/external/decorator/decorator3/_decorator3.py
./networkx/external/decorator/decorator3/__init__.py
./networkx/external/decorator/__init__.py
./networkx/external/decorator/decorator2/_decorator2.py
./networkx/external/decorator/decorator2/__init__.py
./networkx/external/__init__.py
./networkx/algorithms/richclub.py
./networkx/algorithms/assortativity/tests/test_pairs.py
./networkx/algorithms/assortativity/tests/test_correlation.py
./networkx/algorithms/assortativity/tests/test_neighbor_degree.py
./networkx/algorithms/assortativity/tests/test_mixing.py
./networkx/algorithms/assortativity/tests/base_test.py
./networkx/algorithms/assortativity/tests/test_connectivity.py
./networkx/algorithms/assortativity/pairs.py
./networkx/algorithms/assortativity/correlation.py
./networkx/algorithms/assortativity/mixing.py
./networkx/algorithms/assortativity/__init__.py
./networkx/algorithms/assortativity/neighbor_degree.py                                                                                            
./networkx/algorithms/assortativity/connectivity.py                                                                                               
./networkx/algorithms/link_analysis/tests/test_hits.py                                                                                            
./networkx/algorithms/link_analysis/tests/test_pagerank.py                                                                                        
./networkx/algorithms/link_analysis/pagerank_alg.py                                                                                               
./networkx/algorithms/link_analysis/__init__.py                                                                                                   
./networkx/algorithms/link_analysis/hits_alg.py                                                                                                   
./networkx/algorithms/smetric.py                                                                                                                  
./networkx/algorithms/graphical.py                                                                                                                
./networkx/algorithms/shortest_paths/tests/test_generic.py                                                                                        
./networkx/algorithms/shortest_paths/tests/test_weighted.py                                                                                       
./networkx/algorithms/shortest_paths/tests/test_astar.py                                                                                          
./networkx/algorithms/shortest_paths/tests/test_dense_numpy.py                                                                                    
./networkx/algorithms/shortest_paths/tests/test_dense.py                                                                                          
./networkx/algorithms/shortest_paths/tests/test_unweighted.py                                                                                     
./networkx/algorithms/shortest_paths/dense.py                                                                                                     
./networkx/algorithms/shortest_paths/unweighted.py                                                                                                
./networkx/algorithms/shortest_paths/weighted.py                                                                                                  
./networkx/algorithms/shortest_paths/astar.py                                                                                                     
./networkx/algorithms/shortest_paths/__init__.py                                                                                                  
./networkx/algorithms/shortest_paths/generic.py                                                                                                   
./networkx/algorithms/simple_paths.py                                                                                                             
./networkx/algorithms/vitality.py                                                                                                                 
./networkx/algorithms/tests/test_mis.py                                                                                                           
./networkx/algorithms/tests/test_hierarchy.py                                                                                                     
./networkx/algorithms/tests/test_swap.py                                                                                                          
./networkx/algorithms/tests/test_simple_paths.py                                                                                                  
./networkx/algorithms/tests/test_cluster.py                                                                                                       
./networkx/algorithms/tests/test_boundary.py                                                                                                      
./networkx/algorithms/tests/test_distance_regular.py                                                                                              
./networkx/algorithms/tests/test_euler.py                                                                                                         
./networkx/algorithms/tests/test_graphical.py
./networkx/algorithms/tests/test_dag.py
./networkx/algorithms/tests/test_smetric.py
./networkx/algorithms/tests/test_richclub.py
./networkx/algorithms/tests/test_distance_measures.py
./networkx/algorithms/tests/test_matching.py
./networkx/algorithms/tests/test_block.py
./networkx/algorithms/tests/test_cycles.py
./networkx/algorithms/tests/test_clique.py
./networkx/algorithms/tests/test_vitality.py
./networkx/algorithms/tests/test_mst.py
./networkx/algorithms/tests/test_core.py
./networkx/algorithms/approximation/independent_set.py
./networkx/algorithms/approximation/tests/test_independent_set.py
./networkx/algorithms/approximation/tests/test_vertex_cover.py
./networkx/algorithms/approximation/tests/test_ramsey.py
./networkx/algorithms/approximation/tests/test_dominating_set.py
./networkx/algorithms/approximation/tests/test_matching.py
./networkx/algorithms/approximation/tests/test_clique.py
./networkx/algorithms/approximation/dominating_set.py
./networkx/algorithms/approximation/matching.py
./networkx/algorithms/approximation/ramsey.py
./networkx/algorithms/approximation/vertex_cover.py
./networkx/algorithms/approximation/clique.py
./networkx/algorithms/approximation/__init__.py
./networkx/algorithms/chordal/tests/test_chordal.py
./networkx/algorithms/chordal/chordal_alg.py
./networkx/algorithms/chordal/__init__.py
./networkx/algorithms/hierarchy.py
./networkx/algorithms/block.py
./networkx/algorithms/core.py
./networkx/algorithms/distance_regular.py
./networkx/algorithms/components/tests/test_biconnected.py
./networkx/algorithms/components/tests/test_weakly_connected.py
./networkx/algorithms/components/tests/test_attracting.py
./networkx/algorithms/components/tests/test_connected.py
./networkx/algorithms/components/tests/test_strongly_connected.py
./networkx/algorithms/components/strongly_connected.py
./networkx/algorithms/components/weakly_connected.py
./networkx/algorithms/components/attracting.py
./networkx/algorithms/components/biconnected.py
./networkx/algorithms/components/__init__.py
./networkx/algorithms/components/connected.py
./networkx/algorithms/centrality/closeness.py
./networkx/algorithms/centrality/tests/test_load_centrality.py
./networkx/algorithms/centrality/tests/test_betweenness_centrality_subset.py
./networkx/algorithms/centrality/tests/test_communicability.py
./networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality_subset.py
./networkx/algorithms/centrality/tests/test_katz_centrality.py
./networkx/algorithms/centrality/tests/test_current_flow_closeness.py
./networkx/algorithms/centrality/tests/test_eigenvector_centrality.py
./networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality.py
./networkx/algorithms/centrality/tests/test_closeness_centrality.py
./networkx/algorithms/centrality/tests/test_betweenness_centrality.py
./networkx/algorithms/centrality/tests/test_degree_centrality.py
./networkx/algorithms/centrality/katz.py
./networkx/algorithms/centrality/betweenness.py
./networkx/algorithms/centrality/flow_matrix.py
./networkx/algorithms/centrality/load.py
./networkx/algorithms/centrality/current_flow_closeness.py
./networkx/algorithms/centrality/communicability_alg.py
./networkx/algorithms/centrality/eigenvector.py
./networkx/algorithms/centrality/degree_alg.py
./networkx/algorithms/centrality/current_flow_betweenness_subset.py
./networkx/algorithms/centrality/current_flow_betweenness.py
./networkx/algorithms/centrality/betweenness_subset.py
./networkx/algorithms/centrality/__init__.py
./networkx/algorithms/matching.py
./networkx/algorithms/dag.py
./networkx/algorithms/boundary.py
./networkx/algorithms/euler.py
./networkx/algorithms/flow/tests/test_mincost.py
./networkx/algorithms/flow/tests/test_maxflow.py
./networkx/algorithms/flow/tests/test_maxflow_large_graph.py
./networkx/algorithms/flow/mincost.py
./networkx/algorithms/flow/maxflow.py
./networkx/algorithms/flow/__init__.py
./networkx/algorithms/distance_measures.py
./networkx/algorithms/isolate.py
./networkx/algorithms/community/tests/test_kclique.py
./networkx/algorithms/community/kclique.py
./networkx/algorithms/community/__init__.py
./networkx/algorithms/swap.py
./networkx/algorithms/cycles.py
./networkx/algorithms/connectivity/tests/test_cuts.py
./networkx/algorithms/connectivity/tests/test_connectivity.py
./networkx/algorithms/connectivity/__init__.py
./networkx/algorithms/connectivity/connectivity.py
./networkx/algorithms/connectivity/cuts.py
./networkx/algorithms/clique.py
./networkx/algorithms/__init__.py
./networkx/algorithms/mst.py
./networkx/algorithms/bipartite/centrality.py
./networkx/algorithms/bipartite/projection.py
./networkx/algorithms/bipartite/tests/test_basic.py
./networkx/algorithms/bipartite/tests/test_centrality.py
./networkx/algorithms/bipartite/tests/test_cluster.py
./networkx/algorithms/bipartite/tests/test_spectral_bipartivity.py
./networkx/algorithms/bipartite/tests/test_project.py
./networkx/algorithms/bipartite/basic.py
./networkx/algorithms/bipartite/spectral.py
./networkx/algorithms/bipartite/__init__.py
./networkx/algorithms/bipartite/redundancy.py
./networkx/algorithms/bipartite/cluster.py
./networkx/algorithms/traversal/tests/test_bfs.py
./networkx/algorithms/traversal/tests/test_dfs.py
./networkx/algorithms/traversal/depth_first_search.py
./networkx/algorithms/traversal/__init__.py
./networkx/algorithms/traversal/breadth_first_search.py
./networkx/algorithms/mis.py
./networkx/algorithms/operators/unary.py
./networkx/algorithms/operators/tests/test_product.py
./networkx/algorithms/operators/tests/test_binary.py
./networkx/algorithms/operators/tests/test_all.py
./networkx/algorithms/operators/tests/test_unary.py
./networkx/algorithms/operators/all.py
./networkx/algorithms/operators/product.py
./networkx/algorithms/operators/binary.py
./networkx/algorithms/operators/__init__.py
./networkx/algorithms/isomorphism/isomorphvf2.py
./networkx/algorithms/isomorphism/tests/iso_r01_s80.A99
./networkx/algorithms/isomorphism/tests/test_isomorphvf2.py
./networkx/algorithms/isomorphism/tests/iso_r01_s80.B99
./networkx/algorithms/isomorphism/tests/si2_b06_m200.B99
./networkx/algorithms/isomorphism/tests/test_vf2userfunc.py
./networkx/algorithms/isomorphism/tests/si2_b06_m200.A99
./networkx/algorithms/isomorphism/tests/test_isomorphism.py
./networkx/algorithms/isomorphism/matchhelpers.py
./networkx/algorithms/isomorphism/__init__.py
./networkx/algorithms/isomorphism/isomorph.py
./networkx/algorithms/isomorphism/vf2userfunc.py
./networkx/algorithms/cluster.py
./networkx/tests/test_convert.py
./networkx/tests/test_convert_numpy.py
./networkx/tests/test_exceptions.py
./networkx/tests/test_relabel.py
./networkx/tests/test.py
./networkx/tests/__init__.py
./networkx/tests/test_convert_scipy.py
./networkx/tests/benchmark.py
./networkx/generators/tests/test_random_graphs.py
./networkx/generators/tests/test_small.py
./networkx/generators/tests/test_geometric.py
./networkx/generators/tests/test_stochastic.py
./networkx/generators/tests/test_classic.py
./networkx/generators/tests/test_threshold.py
./networkx/generators/tests/test_atlas.py
./networkx/generators/tests/test_line.py
./networkx/generators/tests/test_random_clustered.py
./networkx/generators/tests/test_hybrid.py
./networkx/generators/tests/test_directed.py
./networkx/generators/tests/test_degree_seq.py
./networkx/generators/tests/test_bipartite.py
./networkx/generators/tests/test_intersection.py
./networkx/generators/tests/test_ego.py
./networkx/generators/classic.py
./networkx/generators/social.py
./networkx/generators/intersection.py
./networkx/generators/threshold.py
./networkx/generators/ego.py
./networkx/generators/atlas.py
./networkx/generators/hybrid.py
./networkx/generators/directed.py
./networkx/generators/random_graphs.py
./networkx/generators/bipartite.py
./networkx/generators/line.py
./networkx/generators/__init__.py
./networkx/generators/geometric.py
./networkx/generators/random_clustered.py
./networkx/generators/small.py
./networkx/generators/degree_seq.py
./networkx/generators/stochastic.py
./networkx/exception.py
./networkx/version.py
./networkx/release.py
./networkx/linalg/tests/test_spectrum.py
./networkx/linalg/tests/test_laplacian.py
./networkx/linalg/tests/test_graphmatrix.py
./networkx/linalg/attrmatrix.py
./networkx/linalg/graphmatrix.py
./networkx/linalg/__init__.py
./networkx/linalg/laplacianmatrix.py
./networkx/linalg/spectrum.py
./networkx/readwrite/leda.py
./networkx/readwrite/tests/test_p2g.py
./networkx/readwrite/tests/test_leda.py
./networkx/readwrite/tests/test_shp.py
./networkx/readwrite/tests/test_gpickle.py
./networkx/readwrite/tests/test_yaml.py
./networkx/readwrite/tests/test_graphml.py
./networkx/readwrite/tests/test_gexf.py
./networkx/readwrite/tests/test_sparsegraph6.py
./networkx/readwrite/tests/test_gml.py
./networkx/readwrite/tests/test_adjlist.py
./networkx/readwrite/tests/test_edgelist.py
./networkx/readwrite/tests/test_pajek.py
./networkx/readwrite/multiline_adjlist.py
./networkx/readwrite/json_graph/tests/test_tree.py
./networkx/readwrite/json_graph/tests/test_adjacency.py
./networkx/readwrite/json_graph/tests/test_serialize.py
./networkx/readwrite/json_graph/tests/test_node_link.py
./networkx/readwrite/json_graph/tree.py
./networkx/readwrite/json_graph/serialize.py
./networkx/readwrite/json_graph/node_link.py
./networkx/readwrite/json_graph/__init__.py
./networkx/readwrite/json_graph/adjacency.py
./networkx/readwrite/edgelist.py
./networkx/readwrite/gml.py
./networkx/readwrite/pajek.py
./networkx/readwrite/graphml.py
./networkx/readwrite/gexf.py
./networkx/readwrite/adjlist.py
./networkx/readwrite/gpickle.py
./networkx/readwrite/__init__.py
./networkx/readwrite/sparsegraph6.py
./networkx/readwrite/nx_yaml.py
./networkx/readwrite/p2g.py
./networkx/readwrite/nx_shp.py
./networkx/relabel.py
./networkx/classes/tests/test_multigraph.py
./networkx/classes/tests/test_graph.py
./networkx/classes/tests/test_multidigraph.py
./networkx/classes/tests/test_digraph.py
./networkx/classes/tests/test_graph_historical.py
./networkx/classes/tests/historical_tests.py
./networkx/classes/tests/test_digraph_historical.py
./networkx/classes/tests/test_function.py
./networkx/classes/graph.py
./networkx/classes/multigraph.py
./networkx/classes/function.py
./networkx/classes/digraph.py
./networkx/classes/__init__.py
./networkx/classes/multidigraph.py
./networkx/utils/tests/test_misc.py
./networkx/utils/tests/test_decorators.py
./networkx/utils/tests/test_random_sequence.py
./networkx/utils/tests/test_rcm.py
./networkx/utils/tests/test.txt
./networkx/utils/union_find.py
./networkx/utils/decorators.py
./networkx/utils/rcm.py
./networkx/utils/misc.py
./networkx/utils/random_sequence.py
./networkx/utils/__init__.py
./networkx/drawing/tests/test_agraph.py
./networkx/drawing/tests/test_pydot.py
./networkx/drawing/tests/test_layout.py
./networkx/drawing/tests/test_pylab.py
./networkx/drawing/nx_agraph.py
./networkx/drawing/nx_pydot.py
./networkx/drawing/__init__.py
./networkx/drawing/layout.py
./networkx/drawing/nx_pylab.py
./networkx/__init__.py
./networkx/convert.py
./networkx/testing/tests/test_utils.py
./networkx/testing/utils.py
./networkx/testing/__init__.py

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

./networkx/external/decorator/decorator3/_decorator3.py
./networkx/external/decorator/decorator3/__init__.py
./networkx/external/decorator/__init__.py
./networkx/external/decorator/decorator2/_decorator2.py
./networkx/external/decorator/decorator2/__init__.py
./networkx/external/__init__.py
./networkx/algorithms/richclub.py
./networkx/algorithms/assortativity/tests/test_pairs.py
./networkx/algorithms/assortativity/tests/test_correlation.py
./networkx/algorithms/assortativity/tests/test_neighbor_degree.py
./networkx/algorithms/assortativity/tests/test_mixing.py
./networkx/algorithms/assortativity/tests/base_test.py
./networkx/algorithms/assortativity/tests/test_connectivity.py
./networkx/algorithms/assortativity/pairs.py
./networkx/algorithms/assortativity/correlation.py
./networkx/algorithms/assortativity/mixing.py
./networkx/algorithms/assortativity/__init__.py
./networkx/algorithms/assortativity/neighbor_degree.py                                                                                            
./networkx/algorithms/assortativity/connectivity.py                                                                                               
./networkx/algorithms/link_analysis/tests/test_hits.py                                                                                            
./networkx/algorithms/link_analysis/tests/test_pagerank.py                                                                                        
./networkx/algorithms/link_analysis/pagerank_alg.py                                                                                               
./networkx/algorithms/link_analysis/__init__.py                                                                                                   
./networkx/algorithms/link_analysis/hits_alg.py                                                                                                   
./networkx/algorithms/smetric.py                                                                                                                  
./networkx/algorithms/graphical.py                                                                                                                
./networkx/algorithms/shortest_paths/tests/test_generic.py                                                                                        
./networkx/algorithms/shortest_paths/tests/test_weighted.py                                                                                       
./networkx/algorithms/shortest_paths/tests/test_astar.py                                                                                          
./networkx/algorithms/shortest_paths/tests/test_dense_numpy.py                                                                                    
./networkx/algorithms/shortest_paths/tests/test_dense.py                                                                                          
./networkx/algorithms/shortest_paths/tests/test_unweighted.py                                                                                     
./networkx/algorithms/shortest_paths/dense.py                                                                                                     
./networkx/algorithms/shortest_paths/unweighted.py                                                                                                
./networkx/algorithms/shortest_paths/weighted.py                                                                                                  
./networkx/algorithms/shortest_paths/astar.py                                                                                                     
./networkx/algorithms/shortest_paths/__init__.py                                                                                                  
./networkx/algorithms/shortest_paths/generic.py                                                                                                   
./networkx/algorithms/simple_paths.py                                                                                                             
./networkx/algorithms/vitality.py                                                                                                                 
./networkx/algorithms/tests/test_mis.py                                                                                                           
./networkx/algorithms/tests/test_hierarchy.py                                                                                                     
./networkx/algorithms/tests/test_swap.py                                                                                                          
./networkx/algorithms/tests/test_simple_paths.py                                                                                                  
./networkx/algorithms/tests/test_cluster.py                                                                                                       
./networkx/algorithms/tests/test_boundary.py                                                                                                      
./networkx/algorithms/tests/test_distance_regular.py                                                                                              
./networkx/algorithms/tests/test_euler.py                                                                                                         
./networkx/algorithms/tests/test_graphical.py
./networkx/algorithms/tests/test_dag.py
./networkx/algorithms/tests/test_smetric.py
./networkx/algorithms/tests/test_richclub.py
./networkx/algorithms/tests/test_distance_measures.py
./networkx/algorithms/tests/test_matching.py
./networkx/algorithms/tests/test_block.py
./networkx/algorithms/tests/test_cycles.py
./networkx/algorithms/tests/test_clique.py
./networkx/algorithms/tests/test_vitality.py
./networkx/algorithms/tests/test_mst.py
./networkx/algorithms/tests/test_core.py
./networkx/algorithms/approximation/independent_set.py
./networkx/algorithms/approximation/tests/test_independent_set.py
./networkx/algorithms/approximation/tests/test_vertex_cover.py
./networkx/algorithms/approximation/tests/test_ramsey.py
./networkx/algorithms/approximation/tests/test_dominating_set.py
./networkx/algorithms/approximation/tests/test_matching.py
./networkx/algorithms/approximation/tests/test_clique.py
./networkx/algorithms/approximation/dominating_set.py
./networkx/algorithms/approximation/matching.py
./networkx/algorithms/approximation/ramsey.py
./networkx/algorithms/approximation/vertex_cover.py
./networkx/algorithms/approximation/clique.py
./networkx/algorithms/approximation/__init__.py
./networkx/algorithms/chordal/tests/test_chordal.py
./networkx/algorithms/chordal/chordal_alg.py
./networkx/algorithms/chordal/__init__.py
./networkx/algorithms/hierarchy.py
./networkx/algorithms/block.py
./networkx/algorithms/core.py
./networkx/algorithms/distance_regular.py
./networkx/algorithms/components/tests/test_biconnected.py
./networkx/algorithms/components/tests/test_weakly_connected.py
./networkx/algorithms/components/tests/test_attracting.py
./networkx/algorithms/components/tests/test_connected.py
./networkx/algorithms/components/tests/test_strongly_connected.py
./networkx/algorithms/components/strongly_connected.py
./networkx/algorithms/components/weakly_connected.py
./networkx/algorithms/components/attracting.py
./networkx/algorithms/components/biconnected.py
./networkx/algorithms/components/__init__.py
./networkx/algorithms/components/connected.py
./networkx/algorithms/centrality/closeness.py
./networkx/algorithms/centrality/tests/test_load_centrality.py
./networkx/algorithms/centrality/tests/test_betweenness_centrality_subset.py
./networkx/algorithms/centrality/tests/test_communicability.py
./networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality_subset.py
./networkx/algorithms/centrality/tests/test_katz_centrality.py
./networkx/algorithms/centrality/tests/test_current_flow_closeness.py
./networkx/algorithms/centrality/tests/test_eigenvector_centrality.py
./networkx/algorithms/centrality/tests/test_current_flow_betweenness_centrality.py
./networkx/algorithms/centrality/tests/test_closeness_centrality.py
./networkx/algorithms/centrality/tests/test_betweenness_centrality.py
./networkx/algorithms/centrality/tests/test_degree_centrality.py
./networkx/algorithms/centrality/katz.py
./networkx/algorithms/centrality/betweenness.py
./networkx/algorithms/centrality/flow_matrix.py
./networkx/algorithms/centrality/load.py
./networkx/algorithms/centrality/current_flow_closeness.py
./networkx/algorithms/centrality/communicability_alg.py
./networkx/algorithms/centrality/eigenvector.py
./networkx/algorithms/centrality/degree_alg.py
./networkx/algorithms/centrality/current_flow_betweenness_subset.py
./networkx/algorithms/centrality/current_flow_betweenness.py
./networkx/algorithms/centrality/betweenness_subset.py
./networkx/algorithms/centrality/__init__.py
./networkx/algorithms/matching.py
./networkx/algorithms/dag.py
./networkx/algorithms/boundary.py
./networkx/algorithms/euler.py
./networkx/algorithms/flow/tests/test_mincost.py
./networkx/algorithms/flow/tests/test_maxflow.py
./networkx/algorithms/flow/tests/test_maxflow_large_graph.py
./networkx/algorithms/flow/mincost.py
./networkx/algorithms/flow/maxflow.py
./networkx/algorithms/flow/__init__.py
./networkx/algorithms/distance_measures.py
./networkx/algorithms/isolate.py
./networkx/algorithms/community/tests/test_kclique.py
./networkx/algorithms/community/kclique.py
./networkx/algorithms/community/__init__.py
./networkx/algorithms/swap.py
./networkx/algorithms/cycles.py
./networkx/algorithms/connectivity/tests/test_cuts.py
./networkx/algorithms/connectivity/tests/test_connectivity.py
./networkx/algorithms/connectivity/__init__.py
./networkx/algorithms/connectivity/connectivity.py
./networkx/algorithms/connectivity/cuts.py
./networkx/algorithms/clique.py
./networkx/algorithms/__init__.py
./networkx/algorithms/mst.py
./networkx/algorithms/bipartite/centrality.py
./networkx/algorithms/bipartite/projection.py
./networkx/algorithms/bipartite/tests/test_basic.py
./networkx/algorithms/bipartite/tests/test_centrality.py
./networkx/algorithms/bipartite/tests/test_cluster.py
./networkx/algorithms/bipartite/tests/test_spectral_bipartivity.py
./networkx/algorithms/bipartite/tests/test_project.py
./networkx/algorithms/bipartite/basic.py
./networkx/algorithms/bipartite/spectral.py
./networkx/algorithms/bipartite/__init__.py
./networkx/algorithms/bipartite/redundancy.py
./networkx/algorithms/bipartite/cluster.py
./networkx/algorithms/traversal/tests/test_bfs.py
./networkx/algorithms/traversal/tests/test_dfs.py
./networkx/algorithms/traversal/depth_first_search.py
./networkx/algorithms/traversal/__init__.py
./networkx/algorithms/traversal/breadth_first_search.py
./networkx/algorithms/mis.py
./networkx/algorithms/operators/unary.py
./networkx/algorithms/operators/tests/test_product.py
./networkx/algorithms/operators/tests/test_binary.py
./networkx/algorithms/operators/tests/test_all.py
./networkx/algorithms/operators/tests/test_unary.py
./networkx/algorithms/operators/all.py
./networkx/algorithms/operators/product.py
./networkx/algorithms/operators/binary.py
./networkx/algorithms/operators/__init__.py
./networkx/algorithms/isomorphism/isomorphvf2.py
./networkx/algorithms/isomorphism/tests/iso_r01_s80.A99
./networkx/algorithms/isomorphism/tests/test_isomorphvf2.py
./networkx/algorithms/isomorphism/tests/iso_r01_s80.B99
./networkx/algorithms/isomorphism/tests/si2_b06_m200.B99
./networkx/algorithms/isomorphism/tests/test_vf2userfunc.py
./networkx/algorithms/isomorphism/tests/si2_b06_m200.A99
./networkx/algorithms/isomorphism/tests/test_isomorphism.py
./networkx/algorithms/isomorphism/matchhelpers.py
./networkx/algorithms/isomorphism/__init__.py
./networkx/algorithms/isomorphism/isomorph.py
./networkx/algorithms/isomorphism/vf2userfunc.py
./networkx/algorithms/cluster.py
./networkx/tests/test_convert.py
./networkx/tests/test_convert_numpy.py
./networkx/tests/test_exceptions.py
./networkx/tests/test_relabel.py
./networkx/tests/test.py
./networkx/tests/__init__.py
./networkx/tests/test_convert_scipy.py
./networkx/tests/benchmark.py
./networkx/generators/tests/test_random_graphs.py
./networkx/generators/tests/test_small.py
./networkx/generators/tests/test_geometric.py
./networkx/generators/tests/test_stochastic.py
./networkx/generators/tests/test_classic.py
./networkx/generators/tests/test_threshold.py
./networkx/generators/tests/test_atlas.py
./networkx/generators/tests/test_line.py
./networkx/generators/tests/test_random_clustered.py
./networkx/generators/tests/test_hybrid.py
./networkx/generators/tests/test_directed.py
./networkx/generators/tests/test_degree_seq.py
./networkx/generators/tests/test_bipartite.py
./networkx/generators/tests/test_intersection.py
./networkx/generators/tests/test_ego.py
./networkx/generators/classic.py
./networkx/generators/social.py
./networkx/generators/intersection.py
./networkx/generators/threshold.py
./networkx/generators/ego.py
./networkx/generators/atlas.py
./networkx/generators/hybrid.py
./networkx/generators/directed.py
./networkx/generators/random_graphs.py
./networkx/generators/bipartite.py
./networkx/generators/line.py
./networkx/generators/__init__.py
./networkx/generators/geometric.py
./networkx/generators/random_clustered.py
./networkx/generators/small.py
./networkx/generators/degree_seq.py
./networkx/generators/stochastic.py
./networkx/exception.py
./networkx/version.py
./networkx/release.py
./networkx/linalg/tests/test_spectrum.py
./networkx/linalg/tests/test_laplacian.py
./networkx/linalg/tests/test_graphmatrix.py
./networkx/linalg/attrmatrix.py
./networkx/linalg/graphmatrix.py
./networkx/linalg/__init__.py
./networkx/linalg/laplacianmatrix.py
./networkx/linalg/spectrum.py
./networkx/readwrite/leda.py
./networkx/readwrite/tests/test_p2g.py
./networkx/readwrite/tests/test_leda.py
./networkx/readwrite/tests/test_shp.py
./networkx/readwrite/tests/test_gpickle.py
./networkx/readwrite/tests/test_yaml.py
./networkx/readwrite/tests/test_graphml.py
./networkx/readwrite/tests/test_gexf.py
./networkx/readwrite/tests/test_sparsegraph6.py
./networkx/readwrite/tests/test_gml.py
./networkx/readwrite/tests/test_adjlist.py
./networkx/readwrite/tests/test_edgelist.py
./networkx/readwrite/tests/test_pajek.py
./networkx/readwrite/multiline_adjlist.py
./networkx/readwrite/json_graph/tests/test_tree.py
./networkx/readwrite/json_graph/tests/test_adjacency.py
./networkx/readwrite/json_graph/tests/test_serialize.py
./networkx/readwrite/json_graph/tests/test_node_link.py
./networkx/readwrite/json_graph/tree.py
./networkx/readwrite/json_graph/serialize.py
./networkx/readwrite/json_graph/node_link.py
./networkx/readwrite/json_graph/__init__.py
./networkx/readwrite/json_graph/adjacency.py
./networkx/readwrite/edgelist.py
./networkx/readwrite/gml.py
./networkx/readwrite/pajek.py
./networkx/readwrite/graphml.py
./networkx/readwrite/gexf.py
./networkx/readwrite/adjlist.py
./networkx/readwrite/gpickle.py
./networkx/readwrite/__init__.py
./networkx/readwrite/sparsegraph6.py
./networkx/readwrite/nx_yaml.py
./networkx/readwrite/p2g.py
./networkx/readwrite/nx_shp.py
./networkx/relabel.py
./networkx/classes/tests/test_multigraph.py
./networkx/classes/tests/test_graph.py
./networkx/classes/tests/test_multidigraph.py
./networkx/classes/tests/test_digraph.py
./networkx/classes/tests/test_graph_historical.py
./networkx/classes/tests/historical_tests.py
./networkx/classes/tests/test_digraph_historical.py
./networkx/classes/tests/test_function.py
./networkx/classes/graph.py
./networkx/classes/multigraph.py
./networkx/classes/function.py
./networkx/classes/digraph.py
./networkx/classes/__init__.py
./networkx/classes/multidigraph.py
./networkx/utils/tests/test_misc.py
./networkx/utils/tests/test_decorators.py
./networkx/utils/tests/test_random_sequence.py
./networkx/utils/tests/test_rcm.py
./networkx/utils/tests/test.txt
./networkx/utils/union_find.py
./networkx/utils/decorators.py
./networkx/utils/rcm.py
./networkx/utils/misc.py
./networkx/utils/random_sequence.py
./networkx/utils/__init__.py
./networkx/drawing/tests/test_agraph.py
./networkx/drawing/tests/test_pydot.py
./networkx/drawing/tests/test_layout.py
./networkx/drawing/tests/test_pylab.py
./networkx/drawing/nx_agraph.py
./networkx/drawing/nx_pydot.py
./networkx/drawing/__init__.py
./networkx/drawing/layout.py
./networkx/drawing/nx_pylab.py
./networkx/__init__.py
./networkx/convert.py
./networkx/testing/tests/test_utils.py
./networkx/testing/utils.py
./networkx/testing/__init__.py