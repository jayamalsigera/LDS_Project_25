%% Assert that sensor node i is locally observable, i.e. (A, C_i) is observable.
%
% Local observability requires that each sensor node can individually
% reconstruct the full state from its own output sequence, independently of
% other nodes. The PBH rank test is used: (A, C_i) is observable iff
% rank([lambda*I - A; C_i]) == n for every eigenvalue lambda of A.
%
% In a distributed filter, local observability is a stronger condition than
% collective observability (assertDetectable checks the latter). A node that
% is not locally observable must rely on information from neighbours to
% track all state modes.
%
% Parameters:
%   A      - n x n state transition matrix
%   C_i    - p_i x n output matrix for node i
%   nodeId - (optional) node index, used only in printed messages
%   tol    - (optional) rank tolerance, default 1e-9
%
function assertLocallyObservable(A, C_i, nodeId, tol)
  if nargin < 3 || isempty(nodeId)
    nodeId = [];
  end
  if nargin < 4
    tol = 1e-9;
  end

  n = size(A, 1);
  eigs_A = eig(A);
  unobservable = false;

  for k = 1:numel(eigs_A)
    lambda = eigs_A(k);
    M = [lambda * eye(n) - A; C_i];
    if rank(M, tol) < n
      unobservable = true;
      if isempty(nodeId)
        warning('assertLocallyObservable:notObservable', ...
                'Pair (A, C_i) is NOT locally observable: unobservable mode at lambda = %s.', ...
                num2str(lambda));
      else
        warning('assertLocallyObservable:notObservable', ...
                'Node %d: pair (A, C_%d) is NOT locally observable: unobservable mode at lambda = %s.', ...
                nodeId, nodeId, num2str(lambda));
      end
    end
  end

  if ~unobservable
    if isempty(nodeId)
      fprintf('Pair (A, C_i) is locally observable.\n');
    else
      fprintf('Node %d: pair (A, C_%d) is locally observable.\n', nodeId, nodeId);
    end
  end
end
