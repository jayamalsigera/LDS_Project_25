%% Nominal Distributed Kalman Filter with Stochastic Trigger (SDKF)
%
% The DKF from Battistelli et al. (2018) with the deterministic trigger
% replaced by the stochastic event-triggered schedule from Han et al. (2015).
% No robust (minimax) treatment — b = 0, so Psi = Omega everywhere.
%
% (Han et al., 2015) consider the sensor nodes primitive and the trigger is for
% communicating with the estimator. (Battistely et al., 2018) assumes that each
% node is capable of running the estimation logic, so the trigger is used to
% propagate the state estimates.
%
% This serves as a direct comparison partner for SRDKF: same stochastic
% trigger, same network fusion, but no model-uncertainty robustness.
%
classdef SDKF
  properties
    Ts
    T
    % Tuning parameters
    Z
    % Network
    G
    N
    S
    Pi
    % State estimate history
    X_hat
    % Model matrices
    A
    C
    Q
    R
    n
    senBlock
    % Stats
    RMSE
    txRate
  end

  methods
    function self = SDKF(plant, Ts, T, G, Z)
      self.Ts = Ts;
      self.T  = T;
      self.G  = G;

      self.Z     = Z;

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
        'SDKF:senBlock', 'plant.p (%d) is not divisible by sensor count (%d).', ...
        plant.p, self.S);
      self.senBlock = plant.p / self.S;

      self.X_hat  = zeros(self.n, self.N, T + 1);
      self.txRate = zeros(self.T + 1, 1);
    end

    %% Estimation
    function self = estimate(self, x0_hat, P0, X, Y)
      q_pred     = repmat(P0 \ x0_hat, 1, self.N);
      Omega_pred = repmat(P0 \ eye(self.n), 1, 1, self.N);

      self.X_hat(:, :, 1) = repmat(x0_hat, 1, self.N);

      q_bar     = nan(self.n, self.N);
      Omega_bar = nan(self.n, self.n, self.N);

      for t = 2:self.T + 1
        y = Y(:, t);
        [q_upd, Omega_upd] = self.update(q_pred, Omega_pred, y);

        for i = 1:self.N
          self.X_hat(:, i, t) = pinv(Omega_upd(:, :, i)) * q_upd(:, i);
        end

        c_t = self.exchange(self.X_hat(:, :, t), q_bar, Omega_bar);
        self.txRate(t) = sum(c_t) / self.N;

        [q_fused, Omega_fused] = self.fusion(c_t, q_upd, Omega_upd, q_bar, Omega_bar);
        [q_pred, Omega_pred]   = self.getLocalPriors(q_fused, Omega_fused);
        [q_bar, Omega_bar]     = self.updateGlobalPriors(c_t, q_upd, Omega_upd, q_bar, Omega_bar);
      end

      self.RMSE = calculateRmse(self.X_hat, X);
    end

    %% Correction step — full local measurement update (own measurement always used)
    function [q_upd, Omega_upd] = update(self, q_pred, Omega_pred, y)
      q_upd     = zeros(self.n, self.N);
      Omega_upd = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        if self.G.Nodes(i, :).isSensor
          idx = self.senBlock * (i - 1) + (1:self.senBlock);
          y_i   = y(idx);
          C_i   = self.C(idx, :);
          R_i   = self.R(idx, idx);

          q_upd(:, i)       = q_pred(:, i)       + C_i' * (R_i \ y_i);
          Omega_upd(:, :, i) = Omega_pred(:, :, i) + C_i' * (R_i \ C_i);
        else
          q_upd(:, i)       = q_pred(:, i);
          Omega_upd(:, :, i) = Omega_pred(:, :, i);
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

    %% Prediction step (no robust treatment)
    function [q_pred, Omega_pred] = getLocalPriors(self, q_fused, Omega_fused)
      q_pred     = zeros(self.n, self.N);
      Omega_pred = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        Omega_pred(:, :, i) = predictOmega(Omega_fused(:, :, i), self.A, self.Q);
        q_pred(:, i)        = Omega_pred(:, :, i) * self.A * (Omega_fused(:, :, i) \ q_fused(:, i));
      end
    end

    %% Propagate global priors for non-transmitting nodes
    function [q_bar_next, Omega_bar_next] = updateGlobalPriors(self, c_t, q_upd, Omega_upd, q_bar, Omega_bar)
      q_bar_next     = zeros(self.n, self.N);
      Omega_bar_next = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        if c_t(i)
          q_check     = q_upd(:, i);
          Omega_check = Omega_upd(:, :, i);
        else
          q_check     = q_bar(:, i);
          Omega_check = Omega_bar(:, :, i);
        end

        Omega_bar_next(:, :, i) = predictOmega(Omega_check, self.A, self.Q);
        q_bar_next(:, i)        = Omega_bar_next(:, :, i) * self.A * (Omega_check \ q_check);
      end
    end
  end
end
