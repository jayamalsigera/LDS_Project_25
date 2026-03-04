%% Calculate weights pi for a graph
%
% Inputs:
%   G - graph object (directed or undirected)
%
% Outputs:
%   pi - NxN weight matrix where N is the number of nodes
%
% The Weights pi used in the Ghion, Zorzi (2023) are calculated as:
% pi_ij = (d_i+1)^-1 if i and j are connected
% pi_ij = 0 if i and j are not connected
% pi_ii = (d_i+1)^-1 
% 
% TODO: Review algorithm
function pi = calcPi(G)
    N = numnodes(G);
    pi = zeros(N, N);

    % Get adjacency matrix
    A = adjacency(G, 'weighted');

    % Calculate degree for each node
    degrees = sum(A, 2);

    % Calculate weights for connected nodes
    for i = 1:N
        for j = 1:N
            if i ~= j && A(i, j) > 0   
                pi(i, j) = 1 / (1 + degrees(i));
            end
        end
    end

    %And the self-weights
    for i = 1:N
        pi(i, i) = 1 / (1 + degrees(i));
    end
end
