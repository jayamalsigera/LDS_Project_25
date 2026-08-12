%% Report how many sensor nodes are locally observable, i.e. (A, C_i) observable.
%
% Parameters:
%   A           - n x n state transition matrix
%   C           - p x n stacked output matrix for all sensor nodes
%   sensorCount - number of sensor nodes (C is split into equal row blocks)
%   tol         - (optional) rank tolerance, default 1e-9
%
function assertLocallyObservable(A, C, sensorCount, tol)
  if nargin < 4
    tol = 1e-9;
  end

  n = size(A, 1);
  eigs_A = eig(A);
  senBlock = size(C, 1) / sensorCount;   % measurement rows per sensor
  observableCount = 0;

  for i = 1:sensorCount
    idx = senBlock * (i - 1) + (1:senBlock);
    C_i = C(idx, :);

    observable = true;
    for k = 1:numel(eigs_A)
      M = [eigs_A(k) * eye(n) - A; C_i];
      if rank(M, tol) < n
        observable = false;
        break;
      end
    end

    observableCount = observableCount + observable;
  end

  fprintf('Locally observable nodes: %d/%d (%.0f%%)\n', ...
          observableCount, sensorCount, 100 * observableCount / sensorCount);
end
