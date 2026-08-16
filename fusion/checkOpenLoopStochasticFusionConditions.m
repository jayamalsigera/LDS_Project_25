%% Open-loop stochastic event-triggered transmission criterion.
%
% Based on Han et al. (2015), eq. (8) and (9).
%
% The trigger depends only on the current raw measurement y_k — no feedback
% from the estimator is used. The probability of NOT transmitting is
% mu(y) = exp(-1/2 * y' * Z * y), which is small when y is large (i.e. an
% informative measurement is more likely to be transmitted).
%
% An i.i.d. uniform random variable zeta is drawn; the node transmits iff
% zeta > mu(y), i.e. with probability 1 - mu(y).
%
% Inputs:
%   y  - local raw measurement at sensor node i  (m x 1)
%   Z  - positive-definite weight matrix in measurement space  (m x m)
%
% Output:
%   c  - 1 to transmit, 0 otherwise
%
function c = checkOpenLoopStochasticFusionConditions(y, Z)

  nu = exp(-1/2 * y' * Z * y);

  zeta = rand;

  if zeta > nu
    c = 1;
  else
    c = 0;
  end

end
