%% Calculate Metropolis-Hastings weights for a graph
%
% The Metropolis weights are defined as:
%   W(i,j) = 1 / (1 + max(in_deg(i), out_deg(j))) if j is an IN-neighbor of i
%   W(i,i) = 1 - sum(W(i,j) for j~=i)
%
% IMPORTANT — edge direction. `DKF.fusion`, `SDKF` and `DSEACP` all fuse over
% `inedges(G, i)`, i.e. over the set of nodes that transmit *to* i. The weights
% must be built over that same set or the row sums seen by the fusion loop fall
% below 1 and information leaks away silently on any directed graph. Furthermore,
% on directed graphs we explicitly compare the receiver's in-degree with the
% sender's out-degree.
%
% Note: The degree count includes the self-loop that every node in these graphs
% carries, so it is deg+1 relative to the textbook Metropolis definition. That is
% pre-existing behaviour and is left unchanged here so symmetric-graph results
% stay bit-for-bit identical; only the edge direction handling is corrected.
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

    % On a directed graph we must differentiate between in- and out-degrees.
    % We use A > 0 to count topological edges rather than summing edge weights.
    % Note: this still includes the self-loop if the graph carries one,
    % preserving the pre-existing behaviour for symmetric graphs.
    in_degrees = sum(A > 0, 1)';
    out_degrees = sum(A > 0, 2);

    % Calculate weights for the in-neighbors of each node
    for i = 1:N
        for j = 1:N
            if i ~= j && A(j, i) > 0
                % Metropolis weight for the edge j -> i on a directed graph:
                % we compare the receiver's in-degree with the sender's out-degree.
                max_degree = max(in_degrees(i), out_degrees(j));
                Pi(i, j) = 1 / (1 + max_degree);
            end
        end
    end

    % Set diagonal (self-loop weights) to ensure row sums equal 1
    for i = 1:N
        Pi(i, i) = 1 - sum(Pi(i, 1:end ~= i));
    end
end
