%% Calculate the density of a graph as a percentage.
%
% Returns the percentage of actual edges compared to the maximum possible
% edges (for an undirected graph without self-loops).
%
% Parameters:
%   G - A graph or digraph object
%
% Returns:
%   percentage - Percentage of connections (0-100)
%
function percentage = getConnectivityPercentage(G)
    N = numnodes(G);
    maxPossibleEdges = N * (N - 1);  % Directed graph without self-loops
    percentage = (numedges(G) / maxPossibleEdges) * 100;
end
