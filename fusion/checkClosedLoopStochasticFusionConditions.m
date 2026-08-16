%% Closed-loop stochastic event-triggered transmission criterion.
%
% Based on Han et al. (2015), eq. (8) and (10).
%
% The measurement-space innovation is z_k = y_k - C_i * x_bar, where
% x_bar = (Psi_bar)^{-1} q_bar is the predicted state from the most
% recently transmitted pair (the globally-known a priori prediction).
% This matches Han et al.'s z_k = y_k - y_hat_k^- where y_hat_k^- = C x_hat_k^-.
%
% Using measurement-space keeps Z in S^m_{++} (consistent with the paper)
% and avoids using the posterior state estimate (which would be circular:
% deciding whether to transmit the correction using the correction itself).
%
% Inputs:
%   y_i   - local raw measurement at sensor node i  (m x 1)
%   C_i   - local observation matrix at node i  (m x n)
%   x_bar - predicted state from last transmitted pair  (n x 1)
%   Z     - positive-definite weight matrix in measurement space  (m x m)
%
% Output:
%   c     - 1 to transmit, 0 otherwise
function c = checkClosedLoopStochasticFusionConditions(y_i, C_i, x_bar, Z)

  z = y_i - C_i * x_bar;
  nu = exp(-1/2 * z' * Z * z);

  zeta = rand;

  if zeta > nu
    c = 1;
  else
    c = 0;
  end

end
