% Graph Laplacian eigenvalue optimization
% From a talk delivered by S. Boyd to the 2006 ICM:
% http://www.stanford.edu/~boyd/cvx_opt_graph_lapl_eigs.html
%
%  larger_example   - FDLA and FMMC solutions for a 50-node, 200-edge graph
%  cut_grid_example - FDLA and FMMC solutions for a 64-node, 95-edge cut-grid graph
%  small_example    - FDLA and FMMC solutions for an 8-node, 13-edge graph
%  fmmc             - Computes fastest mixing Markov chain (FMMC) edge weights
%  mh               - Computes the Metropolis-Hastings heuristic edge weights
%  best_const       - Computes the constant edge weight that yields fastest averaging.
%  fdla             - Computes the fastest distributed linear averaging (FDLA) edge weights
%  max_deg          - Computes the maximum-degree heuristic edge weights
%  cut_grid_data    - Generate a cut-grid graph for the ICM 2006 talk example
%  plotgraph        - Plots a graph with each edge width proportional to its weight.
help Contents
