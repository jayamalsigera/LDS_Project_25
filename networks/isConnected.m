%% Check whether a graph is connected.
function C = isConnected(G)
  bins = conncomp(G, 'Type', 'strong');
  nComponents = max(bins);
  C = nComponents == 1;
end
