%% Calculate Fusion Weights (Pi) for a particular graph
%
% The weights are computed according to Ghion, Zorzi (2023):
%   pi_ij = (d_i+1)^-1 if i and j are connected
%   pi_ij = 0 if i and j are not connected
%   pi_ii = (d_i+1)^-1
%
% where d_i is the number of distinct neighbors of i, EXCLUDING the self-loop.
% Each node fuses over its d_i neighbors plus itself (d_i+1 terms), so counting
% self in the degree would leave every row summing to d_i/(d_i+1) < 1 (a
% sub-stochastic, information-leaking consensus). We therefore strip the
% self-loop before computing the degree so rows sum to exactly 1.
%
% This algorithm is very similar to Metropolis-Hastings, but it only looks for
% the degrees of i, not j.
%
% Inputs:
%   G - graph object (directed or undirected)
%
% Outputs:
%   Pi - NxN weight matrix where N is the number of nodes
%
function Pi = calcFusionWeights(G)
    N = numnodes(G);
    Pi = zeros(N, N);

    % Get adjacency matrix
    A = adjacency(G, 'weighted');

    % Calculate degree for each node as the number of distinct neighbors,
    % excluding the self-loop (the graph carries one on every node). Counting
    % it would make each row sub-stochastic; see the header note.
    Aoff = A - diag(diag(A));
    degrees = sum(Aoff > 0, 2);

    % Calculate weights for connected nodes
    for i = 1:N
        for j = 1:N
            if i ~= j && A(i, j) > 0
                Pi(i, j) = 1 / (1 + degrees(i));
            end
        end
    end

    % And the self-weights
    for i = 1:N
        Pi(i, i) = 1 / (1 + degrees(i));
    end
end
