%% Calculate Fusion Weights (Pi) for a particular graph
%
% The weights are computed according to Ghion, Zorzi (2023):
%   pi_ij = (d_i+1)^-1 if j is an IN-neighbor of i (i.e. j in N_i)
%   pi_ij = 0 otherwise
%   pi_ii = (d_i+1)^-1
%
% where d_i is the number of distinct IN-neighbors of i (the nodes i receives
% from), EXCLUDING the self-loop. Each node fuses over its d_i in-neighbors plus
% itself (d_i+1 terms), so counting self in the degree would leave every row
% summing to d_i/(d_i+1) < 1 (a sub-stochastic, information-leaking consensus).
% We therefore strip the self-loop before computing the degree so rows sum to
% exactly 1.
%
% IMPORTANT — edge direction. The estimators fuse over `inedges(G, i)` (see
% `RDKF.fusion`, `DKF.fusion`, `SDKF`, `SRDKF`, `DSEACP`), i.e. over the set
% N_i of nodes that transmit *to* i. The weights must therefore be built from
% the in-neighborhood, matching the paper's Pi (it defines pi_ij for j in N_i).
% On a symmetric graph in- and out-neighborhoods coincide, so the distinction is
% invisible; on a genuine digraph, building Pi over out-neighbors makes the row
% sums over the set `fusion()` actually visits fall well below 1 (measured
% median 0.2 on an ER 4% digraph), and the fused information decays
% geometrically to zero with no error raised. See
% docs/rdkf-tx-saturation-analysis.md §6.4.1.
%
% This algorithm is very similar to Metropolis-Hastings, but it only looks for
% the degrees of i, not j.
%
% Inputs:
%   G - graph object (directed or undirected)
%
% Outputs:
%   Pi - NxN weight matrix where N is the number of nodes. Row i sums to 1 over
%        {i} union N_i (the in-neighborhood), which is exactly the set the
%        estimators' fusion loops iterate.
%
function Pi = calcFusionWeights(G)
    N = numnodes(G);
    Pi = zeros(N, N);

    % Get adjacency matrix. A(u, v) > 0 means there is an edge u -> v, so the
    % IN-neighbors of i are the nonzeros of column i.
    A = adjacency(G, 'weighted');

    % Calculate the in-degree of each node as the number of distinct nodes it
    % receives from, excluding the self-loop (the graph carries one on every
    % node). Counting it would make each row sub-stochastic; see the header.
    Aoff = A - diag(diag(A));
    degrees = sum(Aoff > 0, 1)';   % column sums = in-degree

    % Calculate weights for the in-neighbors of each node
    for i = 1:N
        for j = 1:N
            if i ~= j && A(j, i) > 0
                Pi(i, j) = 1 / (1 + degrees(i));
            end
        end
    end

    % And the self-weights
    for i = 1:N
        Pi(i, i) = 1 / (1 + degrees(i));
    end
end
