%% Check whether a node should transmit according to stochastic triggers
%
% Based on Han et al. (2015), eq. (8) and (9), but using a delta on the state,
% not on the measurements or innovation.
%
% Inputs:
%   stateDelta - state discrepancy used to calculate probability of idle
%   Z          - positive-definite weight matrix in state space (n x n)
%
% Output:
%   c  - 1 to transmit, 0 otherwise
%
function c = checkStochasticTxRule(stateDelta, Z)
  mu = exp(-1/2 * stateDelta' * Z * stateDelta);  % Probability of idle
  c = double(rand > mu);
end
