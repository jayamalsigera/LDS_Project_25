%% Create a random network of Sensors and Extenders.
%
% TODO: There's no iteration limit, so it's possible for the function to hang
% indefinitely.
%
function G = createRandomNetwork(N, S, connTarget, maxLength)
  totalEdges = ceil(connTarget * N * (N - 1));

  diagonal = 1:N+1:N*N;  % Linear Indices of the diagonal elements
  candidateEdges = setdiff(1:N*N, diagonal);  % All edges minus self-loops

  while 1
    A = false(N);  % Initialize empty adjacency matrix
    A(diagonal) = 1;  % Ensure self loops

    sampledIndices = randperm(numel(candidateEdges), totalEdges);
    sampledEdges = candidateEdges(sampledIndices);
    A(sampledEdges) = 1;

    if isConnected(digraph(A))
      break
    end
  end

  % Even though it's not as relevant in this case, still assign coordinates for
  % consistency with plotting functions.
  nodeCoordinates = rand(N, 2) * maxLength;

  % Graph Metadata
  xCoordinates = nodeCoordinates(:, 1);
  yCoordinates = nodeCoordinates(:, 2);

  % Consider the first S nodes as sensors and the rest as regular
  % This makes it easier to match their index with the ones in the SSM output
  isSensor = [true(S, 1); false(N - S, 1)];

  % TODO: A bit of redundancy in the coordinates for now. Need to decide which
  % representation is more convenient outside of the function.
  nodeTable = table(isSensor, xCoordinates, yCoordinates, nodeCoordinates);

  G = digraph(A, nodeTable);
end
