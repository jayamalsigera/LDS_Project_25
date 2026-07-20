function bLocal = computeLocalTolerances(plant, P0, bGlobal, netGraph)
%COMPUTELOCALTOLERANCES  Per-node KL tolerances b^i for RDKFLOC / SRDKFLOC.
%
%   bLocal = computeLocalTolerances(plant, P0, bGlobal, netGraph) returns an
%   N-by-1 vector whose i-th entry is node i's local KL tolerance b^i, computed
%   offline from the GLOBAL least-favorable model of radius bGlobal (Ghion &
%   Zorzi 2023, "A robust strategy with local tolerances", Section 5, eq. 26 in
%   its asymptotic/steady-state form):
%
%     b^i = 1/2 [ logdet(K_i) - logdet(Kbar_i) + tr(Kbar_i / K_i) - (n + p_i) ]
%
%   which is exactly D_KL( N(0, Kbar_i) || N(0, K_i) ). Here K (Kbar) is the
%   asymptotic covariance of the joint variable z_t = (x_{t+1}, y_t) given
%   Y_{t-1} under the pseudo-nominal (least-favorable) density, and the
%   subscript i keeps the n state rows plus node i's own p_i measurement rows,
%   deleting the other nodes' measurement rows.
%
%   The two densities share the same robust filtered error covariance V (the
%   LFM's filtered cov) and differ only in the transition/noise model:
%     - pseudo-nominal: nominal transition, unit-covariance driver;
%     - least-favorable: shifted/scaled driver v = H e + L eps (Kv = L L').
%   Both are produced per step by LeastFavorableModel.precompute, evaluated
%   here at a converged interior index t* (the recursions are stationary in the
%   interior, so a single steady-state slice suffices).
%
%   Slicing convention matches the estimators (RDKF.m:104-107): sensor node i
%   occupies measurement rows (2i-1):(2i) of plant.C, so its slice keeps state
%   rows 1:n and measurement rows n+(2i-1):n+2i. Relay (non-sensor) nodes have
%   no measurement and get the state-only (n) slice (p_i = 0).
%
%   Inputs
%     plant    : SingleTarget2dModel (provides A, C, B, D, n, p, m, T).
%     P0       : prior covariance (seeds the LFM forward sweep).
%     bGlobal  : scalar global KL radius (must be > 0; e.g. lfmKlTolerance).
%     netGraph : network digraph; netGraph.Nodes.isSensor selects the slice.
%
%   Output
%     bLocal   : N-by-1 vector, bLocal(i) = b^i, with 0 < b^i <= bGlobal for
%                sensor nodes (paper: local tolerance never exceeds the global).
%
%   See also LeastFavorableModel, RDKFLOC, SRDKFLOC.

  n = plant.n;
  T = plant.T;
  N = numnodes(netGraph);
  isSensor = netGraph.Nodes.isSensor;

  % Global least-favorable model: precompute() fills V, H, L for t = 0..T.
  lfm = LeastFavorableModel(plant, P0, bGlobal, T);

  A  = lfm.A;   C  = lfm.C;
  Bt = lfm.Btil;  Dt = lfm.Dtil;

  % Steady-state index: forward V converges from t=1, and the backward W->H,L
  % recursion has a boundary at t=T+1, so use a deep interior index and assert
  % both recursions have settled there.
  tStar = round(T / 2);
  assertConverged(lfm.V, tStar,     'V');
  assertConverged(lfm.H, tStar,     'H');
  assertConverged(lfm.L, tStar,     'L');

  V  = lfm.V(:, :, tStar);
  H  = lfm.H(:, :, tStar);
  Kv = lfm.L(:, :, tStar) * lfm.L(:, :, tStar)';   % LF driver covariance

  % Pseudo-nominal joint covariance K of z = (x_{t+1}, y_t): nominal transition
  % driven by unit-covariance noise, state error e ~ N(0, V). (These Px/Py/Pxy
  % blocks are literally what precompute forms at LeastFavorableModel.m:94-96.)
  Px  = symm(A * V * A' + Bt * Bt');
  Py  = symm(C * V * C' + Dt * Dt');
  Pxy = A * V * C' + Bt * Dt';
  K   = [Px, Pxy; Pxy', Py];

  % Least-favorable joint covariance Kbar: transition driven by v = H e + L eps,
  % so x_{t+1} = A xhat + (A + Bt H) e + Bt L eps and y = C xhat + (C + Dt H) e
  % + Dt L eps; take cov over e ~ N(0, V) _|_ eps ~ N(0, I). (H = 0, Kv = I
  % gives Kbar = K, i.e. bGlobal -> 0 gives b^i -> 0.)
  Abar = A + Bt * H;
  Cbar = C + Dt * H;
  Kxx  = symm(Abar * V * Abar' + Bt * Kv * Bt');
  Kyy  = symm(Cbar * V * Cbar' + Dt * Kv * Dt');
  Kxy  = Abar * V * Cbar' + Bt * Kv * Dt';
  Kbar = [Kxx, Kxy; Kxy', Kyy];

  bLocal = zeros(N, 1);
  for i = 1:N
    if isSensor(i)
      pI  = 2;                                   % 2-row measurement block per node
      idx = [1:n, n + (2 * i - 1), n + 2 * i];   % state rows + node i's own rows
    else
      pI  = 0;                                   % relay: state-only slice
      idx = 1:n;
    end

    Ki    = K(idx, idx);
    Kbari = Kbar(idx, idx);

    bLocal(i) = 0.5 * (logdet(Ki) - logdet(Kbari) + ...
                       trace(Ki \ Kbari) - (n + pI));
  end

  % Paper guarantee: 0 < b^i <= bGlobal for every sensor node. Allow a small
  % numerical slack on the upper bound.
  sens = bLocal(isSensor);
  if any(sens <= 0) || any(sens > bGlobal * (1 + 1e-6))
    warning('computeLocalTolerances:OutOfRange', ...
      ['Some sensor b^i fell outside (0, bGlobal]: min=%.3g, max=%.3g, ' ...
       'bGlobal=%.3g. Check steady-state convergence at t*=%d.'], ...
      min(sens), max(sens), bGlobal, tStar);
  end
end

function d = logdet(M)
  % Stable log-determinant of a symmetric positive-definite matrix.
  R = chol(symm(M));
  d = 2 * sum(log(diag(R)));
end

function M = symm(M)
  M = (M + M') / 2;
end

function assertConverged(arr, tStar, name)
  % Confirm the per-step sweep array is stationary around the interior index.
  res = norm(arr(:, :, tStar) - arr(:, :, tStar - 1), 'fro') ...
      + norm(arr(:, :, tStar + 1) - arr(:, :, tStar), 'fro');
  scale = max(norm(arr(:, :, tStar), 'fro'), 1);
  if res / scale > 1e-5
    warning('computeLocalTolerances:NotConverged', ...
      ['%s has not reached steady state at t*=%d (relative step %.3g). ' ...
       'b^i assumes the asymptotic covariances; results may be off.'], ...
      name, tStar, res / scale);
  end
end
