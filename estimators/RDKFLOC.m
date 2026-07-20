%% Robust Distributed Kalman Filter with per-node LOCAL tolerances
%
% Algorithm 2 of Ghion & Zorzi (2023), "A robust strategy with local
% tolerances" (Section 5). Identical to RDKF (Algorithm 1) except each node i
% uses its own KL tolerance b^i in the two robust prediction steps (eqs. 27,
% 28) instead of the shared global tolerance b. The b^i are computed offline
% from the global least-favorable model (see computeLocalTolerances); because
% the local ambiguity set of each node is characterized more accurately, the
% local tolerances are no larger than the global one and typically smaller
% (paper Fig. 4), giving a lower transmission rate at comparable RMSE.
%
% Everything else -- fusion weights, event trigger, correction, the estimate
% loop -- is inherited from RDKF unchanged: RDKFLOC is just RDKF constructed
% with the per-node b^i vector in place of a scalar b.
%
classdef RDKFLOC < RDKF
  methods
    function self = RDKFLOC(plant, Ts, T, G, alpha, beta, delta, bGlobal, P0)
      % bGlobal : scalar global KL radius (lfmKlTolerance, e.g. 0.05).
      % P0      : prior covariance, needed to build the global LFM.
      % The per-node vector is computed here and handed to the RDKF superclass;
      % it lives in the inherited self.b property, so it is inspectable after
      % construction (filter.b).
      bLocal = computeLocalTolerances(plant, P0, bGlobal, G);
      self@RDKF(plant, Ts, T, G, alpha, beta, delta, bLocal);
    end
  end
end
