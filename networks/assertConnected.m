%% Assert that a graph is connected.
%
% Errors if the graph has more than one connected component. A disconnected
% sensor network breaks the consensus-based filters, which assume every node
% can influence every other via some multi-hop path.
%
% Parameters:
%   G - A graph or digraph object
%
function assertConnected(G)
  bins = conncomp(G, 'Type', 'weak');
  nComponents = max(bins);
  if nComponents > 1
    error('assertConnected:disconnected', ...
          'Network is not connected: %d components found.', nComponents);
  end
  fprintf('Network is connected (%d nodes, %d edges).\n', numnodes(G), numedges(G));
end
