function c = checkDeterministicFusionConditions(X, x_bar, Omega, epsilon)
% Deterministic event-triggered transmission criterion (error-norm only).
% Based on Ghion & Zorzi (2023), partial condition from eq. (9).
%
% A node does NOT transmit if the discrepancy between its current state
% estimate X and the globally-known predicted state x_bar is small in the
% Omega-weighted norm:
%
%   ||X - x_bar||^2_Omega = (X - x_bar)' * Omega * (X - x_bar) < epsilon
%
% This is the error-norm half of eq. (9). The full RDKF trigger also
% includes a Loewner sandwich condition on the covariance matrices
% (implemented inline in the RDKF and DKF classes via loewnerBetweenEig).
%
% Inputs:
%   X       - current state estimate  (n x 1)
%   x_bar   - predicted state from last transmitted pair  (n x 1)
%   Omega   - information matrix used as norm weight  (n x n)
%   epsilon - error-norm threshold (alpha in the paper)
%
% Output:
%   c       - 1 to transmit, 0 otherwise

  e = X - x_bar;

  if e' * Omega * e < epsilon
    c = 0;
  else
    c = 1;
  end

end
