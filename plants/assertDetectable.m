%% Assert that a discrete-time pair (A, C) is detectable.
%
% Detectability is the dual of stabilizability: (A, C) is detectable iff
% rank([lambda*I - A; C]) == n for every eigenvalue lambda of A with
% |lambda| >= 1. If any such eigenvalue fails the rank test, that mode is
% unobservable and unstable — the Kalman filter error covariance will not
% converge.
%
% Parameters:
%   A   - nxn state transition matrix
%   C   - pxn output matrix
%   tol - (optional) rank tolerance, default 1e-9
%
function assertDetectable(A, C, tol)
  if nargin < 3
    tol = 1e-9;
  end

  n = size(A, 1);
  eigs_A = eig(A);
  unstable = eigs_A(abs(eigs_A) >= 1 - tol);

  for k = 1:numel(unstable)
    lambda = unstable(k);
    M = [lambda * eye(n) - A; C];
    if rank(M, tol) < n
      error('assertDetectable:undetectable', ...
            'Pair (A, C) is not detectable: unobservable mode at lambda = %s.', ...
            num2str(lambda));
    end
  end

  fprintf('Pair (A, C) is detectable (%d unstable/marginal mode(s)).\n', numel(unstable));
end
