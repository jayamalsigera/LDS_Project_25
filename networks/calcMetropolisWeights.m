%% Calculate Metropolis-Hastings weights for a graph
%
% Inputs:
%   G - graph object (directed or undirected)
%
% Outputs:
%   W - NxN weight matrix where N is the number of nodes
%
% The Metropolis weights are defined as:
%   W(i,j) = 1 / (1 + max(deg(i), deg(j))) if i~j and edge exists
%   W(i,i) = 1 - sum(W(i,j) for j~=i)
% TODO: Review algorithm
function W = calcMetropolisWeights(G)
    N = numnodes(G);
    W = zeros(N, N);

    % Get adjacency matrix
    A = adjacency(G, 'weighted');

    % Calculate degree for each node
    degrees = sum(A, 2);

    % Calculate weights for connected nodes
    for i = 1:N
        for j = 1:N
            if i ~= j && A(i, j) > 0
                % Metropolis weight for edge (i,j)
                max_degree = max(degrees(i), degrees(j));
                W(i, j) = 1 / (1 + max_degree);
            end
        end
    end

    % Set diagonal (self-loop weights) to ensure row sums equal 1
    for i = 1:N
        W(i, i) = 1 - sum(W(i, 1:end ~= i));
    end
end
