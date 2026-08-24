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
    senBlock
    % Stats
    RMSE
    txRate
    theta_hist
    theta_bar_hist
  end

  methods
    function self = SRDKF(plant, Ts, T, G, b, Z)
      self.Ts = Ts;
      self.T = T;
      self.G = G;

      self.N = numnodes(G);
      self.S = sum(G.Nodes.isSensor);

      self.Pi = calcFusionWeights(G);

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      % Measurement rows per sensor (p_i); derived from C. See RDKF.
      assert(mod(plant.p, self.S) == 0, ...
        'SRDKF:senBlock', 'plant.p (%d) is not divisible by sensor count (%d).', ...
        plant.p, self.S);
      self.senBlock = plant.p / self.S;

      % b may be a scalar (uniform tolerance) or an N-vector (per-node local
      % tolerances b^i, for SRDKFLOC). Store as an N-vector so the prediction
      % steps can index self.b(i); a scalar broadcasts to a constant vector,
      % leaving scalar callers unchanged (b = 0 still recovers SDKF per node).
      if isscalar(b), b = repmat(b, self.N, 1); end
      self.b = b;
      self.Z = Z;

      self.X_hat = zeros(self.n, self.N, T + 1);
      self.txRate = zeros(self.T + 1, 1);
      self.theta_hist = zeros(self.N, T + 1);
      self.theta_bar_hist = zeros(self.N, T + 1);
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

        c_t = self.exchange(self.X_hat(:, :, t), q_bar, Psi_bar);
        self.txRate(t) = sum(c_t) / self.N;

        [q_fused, Omega_fused] = self.fusion(c_t, q_upd, Omega_upd, q_bar, Psi_bar);
        [q_pred, Psi, theta_t] = self.getLocalPriors(q_fused, Omega_fused);
        self.theta_hist(:, t) = theta_t;
        [q_bar, Psi_bar, theta_bar_t] = self.updateGlobalPriors(c_t, q_upd, Omega_upd, q_bar, Psi_bar);
        self.theta_bar_hist(:, t) = theta_bar_t;
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
          % C, R, and y are laid out with a senBlock-row block per node index
          % (p_i = senBlock); non-sensor nodes still occupy their block as placeholders.
          idx = self.senBlock * (i - 1) + (1:self.senBlock);

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

    %% Information exchange
    function c_t = exchange(self, X_hat, q_bar, Omega_bar)
      c_t = ones(self.N, 1);
      if any(isnan(Omega_bar), 'all')
        return % Every node transmits on the first iteration
      end

      for i = 1:self.N
        x_bar_i = Omega_bar(:, :, i) \ q_bar(:, i);
        stateDelta = X_hat(:, i) - x_bar_i;
        c_t(i) = checkStochasticTxRule(stateDelta, self.Z);
      end
    end

    %% Information Pair Fusion
    function [q_fused, Omega_fused] = fusion(self, c_t, q_upd, Omega_upd, q_bar, Psi_bar)
      q_fused = zeros(self.n, self.N);
      Omega_fused = zeros(self.n, self.n, self.N);

      Zinv = self.Z \ eye(self.n);

      for i = 1:self.N
        [~, nids] = inedges(self.G, i);

        for j = nids'
          pi_ij = self.Pi(i, j);

          if (i == j) || c_t(j)
            % Node received active transmission, or it's the node's local update
            q_fused(:, i) = q_fused(:, i) + pi_ij * q_upd(:, j);
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + pi_ij * Omega_upd(:, :, j);
          else
            % Unified reconstruction for ANY silent neighbor (sensor or relay):
            % Reconstruct via state-covariance inflation: P_tilde = P_bar_j + Zinv
            P_bar_j     = Psi_bar(:, :, j) \ eye(self.n);
            P_tilde     = P_bar_j + Zinv;
            Omega_tilde = P_tilde \ eye(self.n);

            x_bar_j     = Psi_bar(:, :, j) \ q_bar(:, j);
            q_tilde     = Omega_tilde * x_bar_j;

            q_fused(:, i) = q_fused(:, i) + pi_ij * q_tilde;
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + pi_ij * Omega_tilde;
          end
        end
      end
    end


    %% Prediction Step
    function [q_pred, Psi, theta_vec] = getLocalPriors(self, q_fused, Omega_fused)
      q_pred = zeros(self.n, self.N);
      Psi = zeros(self.n, self.n, self.N);
      theta_vec = zeros(self.N, 1);

      for i = 1:self.N
        q_i_F = q_fused(:, i);
        Omega_i_F = Omega_fused(:, :, i);

        [q_i_pred, Psi_i_pred, ~, theta_i] = predictRobustFusion( ...
          q_i_F, Omega_i_F, self.A, self.Q, self.b(i));

        q_pred(:, i) = q_i_pred;
        Psi(:, :, i) = Psi_i_pred;
        theta_vec(i) = theta_i;
      end
    end

    % Update the global priors used in the fusion rule when there is no transmission
    function [q_bar_next, Psi_bar_next, theta_bar_vec] = updateGlobalPriors(self, c_t, q_upd, Omega_upd, q_bar, Psi_bar)
      q_bar_next = zeros(self.n, self.N);
      Psi_bar_next = zeros(self.n, self.n, self.N);
      theta_bar_vec = zeros(self.N, 1);

      for i = 1:self.N
        if c_t(i)
          q_i_check = q_upd(:, i);
          Omega_i_check = Omega_upd(:, :, i);
        else
          q_i_check = q_bar(:, i);
          Omega_i_check = Psi_bar(:, :, i);
        end

        [q_i_bar, Psi_i_bar, ~, theta_i_bar] = predictRobustFusion( ...
          q_i_check, Omega_i_check, self.A, self.Q, self.b(i));

        q_bar_next(:, i) = q_i_bar;
        Psi_bar_next(:, :, i) = Psi_i_bar;
        theta_bar_vec(i) = theta_i_bar;
      end
    end
  end
end
