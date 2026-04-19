%% Assert that a discrete-time pair (A, B) is stabilizable.
%
% Uses the PBH test: (A, B) is stabilizable iff rank([lambda*I - A, B]) == n
% for every eigenvalue lambda of A with |lambda| >= 1 (unstable or marginal).
% If any such eigenvalue fails the rank test, that mode is uncontrollable
% and unstable - the Kalman filters built on this model will diverge.
%
% Parameters:
%   A   - nxn state transition matrix
%   B   - nxm input matrix
%   tol - (optional) rank tolerance, default 1e-9
%
function assertStabilizable(A, B, tol)
  if nargin < 3
    tol = 1e-9;
  end

  n = size(A, 1);
  eigs_A = eig(A);
  unstable = eigs_A(abs(eigs_A) >= 1 - tol);

  for k = 1:numel(unstable)
    lambda = unstable(k);
    M = [lambda * eye(n) - A, B];
    if rank(M, tol) < n
      error('assertStabilizable:unstabilizable', ...
            'Pair (A, B) is not stabilizable: uncontrollable mode at lambda = %s.', ...
            num2str(lambda));
    end
  end

  fprintf('Pair (A, B) is stabilizable (%d unstable/marginal mode(s)).\n', numel(unstable));
end
