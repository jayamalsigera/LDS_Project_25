%% RMSE across nodes, per time step.
%
% Given per-node estimates X_hat (n x N x T+1) and the true state
% X (n x T+1), returns a (T+1) x 1 vector where entry t is the
% root-mean-square of per-node full-state error norms:
%
%   rmse(t) = sqrt( mean_i ||X_hat(:, i, t) - X(:, t)||^2 )
%
% For a centralized estimator with a single "node", pass X_hat as
% (n x 1 x T+1).
%
function rmse = calculateRmse(X_hat, X)
  T = size(X, 2) - 1;
  rmse = zeros(T + 1, 1);
  for t = 1:T + 1
    err = X_hat(:, :, t) - X(:, t); % implicit expansion over N
    rmse(t) = sqrt(mean(sum(err .^ 2, 1)));
  end
end
