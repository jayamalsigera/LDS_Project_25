%% Create a random spatial network of Sensors and Extenders.
%
% Algorithm:
% - Nodes are placed in a 2D space uniformly at random
% - Create clusters by connecting nearby nodes up to a distance threshold
% - Interactively connect the clusters by finding the closest nodes available
% - Repeat until graph is fully connected.
%
% The graph is constructed as undirected to make it easier to get strong
% connectivity and then it's converted to a digraph.
%
% All self loops are removed.
%
% Parameters:
% - N - Total number of Nodes (Sensors + Extenders)
% - S - Number of Sensors
% - maxLength - Maximum coordinate value (grid is [0, maxLength]×[0, maxLength])
%
% Returns:
% - G - A `digraph` object representing the Network. Use `G.Nodes` to access
%  node metadata, like their coordinates and whether they're a sensor or not.
%
function G = createSpatialNetwork(N, S, maxLength)
  % Place nodes
  nodeCoordinates = rand(N, 2) * maxLength;

  % Compute pairwise distances between all nodes
  distances = pdist2(nodeCoordinates, nodeCoordinates);

  % Connect nearby nodes, up to a threshold
  connectionRadius = 0.12 * maxLength;
  A = sparse(distances < connectionRadius);

  A(1:N + 1:end) = 0; % Remove self-loops

  % Iteratively connect isolated nodes/small clusters until fully connected
  G = graph(A);
  bins = conncomp(G);
  numBins = max(bins);
  while numBins > 1 % While Graph is NOT fully connected
    selectedBin = max(bins);

    % Isolate the bin nodes from the rest of the nodes
    binNodes = find(bins == selectedBin);
    otherNodes = find(bins ~= selectedBin);

    distancesToOtherNodes = distances(binNodes, otherNodes);
    [minDist, linearIdx] = min(distancesToOtherNodes(:));

    [rowIdx, colIdx] = ind2sub(size(distancesToOtherNodes), linearIdx);

    closestBinNode = binNodes(rowIdx);
    closestOtherNode = otherNodes(colIdx);

    G = addedge(G, closestBinNode, closestOtherNode);

    bins = conncomp(G);
    numBins = max(bins);
  end

  % Graph Metadata
  xCoordinates = nodeCoordinates(:, 1);
  yCoordinates = nodeCoordinates(:, 2);

  % Consider the first S nodes as sensors and the rest as regular
  % This makes it easier to match their index with the ones in the SSM output
  isSensor = [true(S, 1); false(N - S, 1)];

  % TODO: A bit of redundancy in the coordinates for now. Need to decide which
  % representation is more convenient outside of the function.
  nodeTable = table(isSensor, xCoordinates, yCoordinates, nodeCoordinates);

  G = digraph(G.adjacency, nodeTable);
end
