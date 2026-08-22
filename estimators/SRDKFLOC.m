%% Robust Distributed KF with Stochastic Trigger and per-node LOCAL tolerances
%
% The local-tolerance counterpart of SRDKF: the stochastic-trigger robust
% distributed filter (Ghion & Zorzi 2023 robust layer + Han et al. 2015
% stochastic trigger), but with each node i using its own KL tolerance b^i
% (Algorithm 2, eqs. 27/28) instead of the shared global b. The b^i are
% computed offline from the global least-favorable model (see
% computeLocalTolerances).
%
% Since SRDKF's transmission is essentially insensitive to b, the payoff here
% is on estimation quality: the smaller, geometry-matched per-node tolerances
% avoid the over-hedging of a uniform global b. Everything except the tolerance
% is inherited from SRDKF unchanged.
%
classdef SRDKFLOC < SRDKF
  methods
    function self = SRDKFLOC(plant, Ts, T, G, bGlobal, Z, P0)
      % bGlobal     : scalar global KL radius (lfmKlTolerance, e.g. 0.05).
      % Z           : stochastic-trigger error-norm weight (as in SRDKF).
      % P0          : prior covariance, needed to build the global LFM. Placed
      %               last so the signature is an additive tail on SRDKF's.
      bLocal = computeLocalTolerances(plant, P0, bGlobal, G);
      self@SRDKF(plant, Ts, T, G, bLocal, Z);
    end
  end
end
