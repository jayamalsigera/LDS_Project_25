%% RMSE across nodes, per time step.
%
% Given stacked per-node estimates X_hat (n x N x T+1) and the true state
% X (n x T+1), returns a (T+1) x 1 vector where entry t is the RMSE across
% all nodes at time t-1.
%
function rmse = calculateRSME(X_hat, X)
  T = size(X, 2) - 1;
  rmse = zeros(T + 1, 1);
  for t = 1:T
    err = X_hat(:, :, t) - X(:, t); % implicit expansion over N
    rmse(t) = sqrt(mean(sum(err .^ 2, 1)));
  end
end
