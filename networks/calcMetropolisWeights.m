%% Calculate Metropolis-Hastings weights for a graph
%
% The Metropolis weights are defined as:
%   W(i,j) = 1 / (1 + max(deg(i), deg(j))) if j is an IN-neighbor of i
%   W(i,i) = 1 - sum(W(i,j) for j~=i)
%
% IMPORTANT — edge direction. `DKF.fusion`, `SDKF` and `DSEACP` all fuse over
% `inedges(G, i)`, i.e. over the set of nodes that transmit *to* i. The weights
% must be built over that same set or the row sums seen by the fusion loop fall
% below 1 and information leaks away silently on any directed graph. In- and
% out-neighborhoods coincide on a symmetric graph, which is why this stayed
% latent — every graph the repo builds is symmetric. Same defect, and same fix,
% as `calcFusionWeights`; see docs/rdkf-tx-saturation-analysis.md §6.4.1, §6.8.
%
% Note: `degrees` counts the self-loop that every node in these graphs carries,
% so it is deg(i)+1 relative to the textbook Metropolis definition. That is
% pre-existing behaviour and is left unchanged here so symmetric-graph results
% stay bit-for-bit identical; only the edge direction is corrected.
%
% Inputs:
%   G - graph object (directed or undirected)
%
% Outputs:
%   Pi - NxN weight matrix where N is the number of nodes. Row i sums to 1 over
%        {i} union N_i (the in-neighborhood), the set the estimators iterate.
%
function Pi = calcMetropolisWeights(G)
    N = numnodes(G);
    Pi = zeros(N, N);

    % Get adjacency matrix. A(u, v) > 0 means an edge u -> v, so the
    % IN-neighbors of i are the nonzeros of column i.
    A = adjacency(G, 'weighted');

    % In-degree of each node (includes the self-loop; see the header note)
    degrees = sum(A, 1)';

    % Calculate weights for the in-neighbors of each node
    for i = 1:N
        for j = 1:N
            if i ~= j && A(j, i) > 0
                % Metropolis weight for the edge j -> i
                max_degree = max(degrees(i), degrees(j));
                Pi(i, j) = 1 / (1 + max_degree);
            end
        end
    end

    % Set diagonal (self-loop weights) to ensure row sums equal 1
    for i = 1:N
        Pi(i, i) = 1 - sum(Pi(i, 1:end ~= i));
    end
end
