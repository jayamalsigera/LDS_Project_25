%% Robust Distributed Kalman Filter With Stochastic Trigger
%
% Implementation based on Ghion, Zorzi (2023) but with stochastic event
% trigger from Han et al. (2015)
%
% Negative-information fusion (Han et al. 2015, Thm 2, eqs 24-26): a node
% always uses its OWN measurement fully. When a sensor neighbor j stays
% silent, the receiving node reconstructs j's contribution with the
% enlarged-noise pseudo-update on j's globally known (robust) prior — R_eff =
% R_j + Z^-1, mean left at the prior, information still tightened — instead of
% the discounted stale prior q_bar/(1+delta). The robust prediction layer
% (predictRobustFusion / predictNoTransmit, KL tolerance b) is unchanged, so
% b = 0 still recovers SDKF. Validated standalone in tests/clsetKfUnitTest.m.
%
classdef SRDKF
  properties
    Ts
    T
    % Tunning Parameters
    alpha
    beta
    delta
    b
    Z
    % Network Graph and Parameters
    G
    Pi
    N
    S
    % State estimate history
    X_hat
    % Model Matrices
    A
    C
    Q
    R
    n
    % Stats
    RMSE
    txRate
    % flag for open_loop ('open') vs closed_loop ('closed')
    triggerMode
  end

  methods
    function self = SRDKF(plant, Ts, T, G, alpha, beta, delta, b, Z, triggerMode)
      self.Ts = Ts;
      self.T = T;
      self.G = G;

      self.alpha = alpha;
      self.beta = beta;
      self.delta = delta;

      self.N = numnodes(G);
      self.S = sum(G.Nodes.isSensor);

      self.Pi = calcFusionWeights(G);

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      % b may be a scalar (uniform tolerance) or an N-vector (per-node local
      % tolerances b^i, for SRDKFLOC). Store as an N-vector so the prediction
      % steps can index self.b(i); a scalar broadcasts to a constant vector,
      % leaving scalar callers unchanged (b = 0 still recovers SDKF per node).
      if isscalar(b), b = repmat(b, self.N, 1); end
      self.b = b;
      self.Z = Z;

      self.X_hat = zeros(self.n, self.N, T + 1);
      self.txRate = zeros(self.T + 1, 1);

      self.triggerMode = triggerMode;
    end

    %% Estimation Method
    function self = estimate(self, x0_hat, P0, X, Y)
      % It's unrealistic to assume all nodes share same initial conditions
      % (Battistelli & Chisci, 2014), specially with "perfect knowledge", but
      % this allow us to get results in similar shape to (Ghion & Zorzi, 2023)
      q_pred = repmat(P0 \ x0_hat, 1, self.N);
      Psi = repmat(P0 \ eye(self.n), 1, 1, self.N);

      self.X_hat(:, :, 1) = repmat(x0_hat, 1, self.N);

      % Initializing the "global" predictions, assuming c_t = 1 for all nodes
      % in the first iteration (i.e. the first fusion step only relies on the
      % local filtered values).
      q_bar = nan(self.n, self.N);
      Psi_bar = nan(self.n, self.n, self.N);

      for t = 2:self.T + 1
        y = Y(:, t);
        [q_upd, Omega_upd] = self.update(q_pred, Psi, y);

        for i = 1:self.N
          % Return to state representation from information form
          self.X_hat(:, i, t) = pinv(Omega_upd(:, :, i)) * q_upd(:, i);
        end

        c_t = self.exchange(q_bar, Psi_bar, y);
        self.txRate(t) = sum(c_t) / self.N;

        [q_fused, Omega_fused] = self.fusion(c_t, q_upd, Omega_upd, q_bar, Psi_bar);
        [q_pred, Psi] = self.getLocalPriors(q_fused, Omega_fused);
        [q_bar, Psi_bar] = self.updateGlobalPriors(c_t, q_upd, Omega_upd, q_bar, Psi_bar);
      end

      self.RMSE = calculateRmse(self.X_hat, X);
    end

    %% Correction/Update/Measurement step
    % Update the local information pair to obtain q_k|k and Omega_k|k of each
    % node
    function [q_upd, Omega_upd] = update(self, q_pred, Psi, y)
      q_upd = zeros(self.n, self.N);
      Omega_upd = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        if self.G.Nodes(i, :).isSensor
          idx = (2 * i - 1):(2 * i);

          y_i = y(idx);
          C_i = self.C(idx, :);
          R_i = self.R(idx, idx); % R (block) diagonal since D (block) diagonal

          q_upd(:, i) = q_pred(:, i) + C_i' * (R_i \ y_i);
          Omega_upd(:, :, i) = Psi(:, :, i) + C_i' * (R_i \ C_i);
        else
          q_upd(:, i) = q_pred(:, i);
          Omega_upd(:, :, i) = Psi(:, :, i);
        end
      end
    end

    %% Information Exchange
    % While we don't have to transmit anything, in this step we're calculating
    % c^i_t for all nodes.
    function c_t = exchange(self, q_bar, Psi_bar, y)
      c_t = ones(self.N, 1);
      if any(isnan(Psi_bar), "all")
        return % Every node transmits on the first iteration
      end

      for i = 1:self.N
        if ~self.G.Nodes(i, :).isSensor
          % Non-sensor nodes have no measurement, so neither stochastic
          % trigger applies. They always share their global prior (c=0).
          c_t(i) = 0;
          continue
        end

        idx = (2 * i - 1):(2 * i);
        y_i = y(idx);

        switch lower(self.triggerMode)
          case 'closed'
            x_bar_i = Psi_bar(:, :, i) \ q_bar(:, i);
            C_i = self.C(idx, :);
            c_t(i) = checkClosedLoopStochasticFusionConditions(y_i, C_i, x_bar_i, self.Z);

          case 'open'
            c_t(i) = checkOpenLoopStochasticFusionConditions(y_i, self.Z);

          otherwise
            error('Unknown triggerMode. Use ''closed'' or ''open''.');
        end
      end
    end

    %% Information Pair Fusion
    function [q_fused, Omega_fused] = fusion(self, c_t, q_upd, Omega_upd, q_bar, Psi_bar)
      q_fused = zeros(self.n, self.N);
      Omega_fused = zeros(self.n, self.n, self.N);
      Zinv = self.Z \ eye(size(self.Z, 1));

      for i = 1:self.N
        [~, nids] = inedges(self.G, i);

        for j = nids'
          pi_ij = self.Pi(i, j);

          if (i == j) || c_t(j)
            % Node i has received from j or this is a self-loop (node has
            % access to its own local info)
            q_fused(:, i) = q_fused(:, i) + pi_ij * q_upd(:, j);
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + pi_ij * Omega_upd(:, :, j);
          elseif self.G.Nodes(j, :).isSensor
            % Silent sensor neighbor: reconstruct its contribution with Han's
            % negative-information update on j's globally known robust prior.
            % Silence => innovation was sub-threshold => mean stays at the
            % prior x_bar_j while the information still tightens by
            % C_j' inv(R_j + inv(Z)) C_j.
            jdx  = (2 * j - 1):(2 * j);
            C_j  = self.C(jdx, :);
            R_j  = self.R(jdx, jdx);
            Reff = R_j + Zinv;
            info = C_j' * (Reff \ C_j);

            x_bar_j     = Psi_bar(:, :, j) \ q_bar(:, j);
            q_recon     = q_bar(:, j)        + info * x_bar_j;
            Omega_recon = Psi_bar(:, :, j)   + info;

            q_fused(:, i) = q_fused(:, i) + pi_ij * q_recon;
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + pi_ij * Omega_recon;
          else
            % Silent non-sensor (relay) neighbor: no measurement, so silence
            % carries no negative information — keep the discounted prior.
            q_tilde = (1 / (1 + self.delta)) * q_bar(:, j);
            Omega_tilde = (1 / (1 + self.delta)) * Psi_bar(:, :, j);

            q_fused(:, i) = q_fused(:, i) + pi_ij * q_tilde;
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + pi_ij * Omega_tilde;
          end
        end
      end
    end

    %% Prediction Step
    function [q_pred, Psi] = getLocalPriors(self, q_fused, Omega_fused)
      q_pred = zeros(self.n, self.N);
      Psi = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        q_i_F = q_fused(:, i);
        Omega_i_F = Omega_fused(:, :, i);

        [q_i_pred, Psi_i_pred, ~, ~] = predictRobustFusion( ...
          q_i_F, Omega_i_F, self.A, self.Q, self.b(i));

        q_pred(:, i) = q_i_pred;
        Psi(:, :, i) = Psi_i_pred;
      end
    end

    % Update the global priors used in the fusion rule when there is no transmission
    function [q_bar_next, Psi_bar_next] = updateGlobalPriors(self, c_t, q_upd, Omega_upd, q_bar, Psi_bar)
      q_bar_next = zeros(self.n, self.N);
      Psi_bar_next = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        if c_t(i)
          q_i_check = q_upd(:, i);
          Omega_i_check = Omega_upd(:, :, i);
        else
          q_i_check = q_bar(:, i);
          Omega_i_check = Psi_bar(:, :, i);
        end

        [q_i_bar, Psi_i_bar, ~, ~] = predictNoTransmit( ...
          q_i_check, Omega_i_check, self.A, self.Q, self.b(i));

        q_bar_next(:, i) = q_i_bar;
        Psi_bar_next(:, :, i) = Psi_i_bar;
      end
    end

    %% Plotting
    function plotTrajectory(self, X)
      % TODO: Would be cool if we could plot P(t) somehow
      % TODO: Restrict axis to ranges of X

      meanX_hat = squeeze(mean(self.X_hat, 2));

      figure
      plot(meanX_hat(3, :), meanX_hat(4, :));
      hold on
      plot(X(3, :), X(4, :));
      hold off
      title(sprintf("SRDKF Estimated Trajectory (%s-loop)", self.triggerMode))
      xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
      ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
      legend({sprintf('SRDKF-%s', self.triggerMode), "Actual Model"})
      grid()
    end
  end
end
