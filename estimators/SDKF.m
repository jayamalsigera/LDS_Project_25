%% Nominal Distributed Kalman Filter with Stochastic Trigger (SDKF)
%
% The DKF from Battistelli et al. (2018) with the deterministic trigger
% replaced by the stochastic event-triggered schedule from Han et al. (2015).
% No robust (minimax) treatment — b = 0, so Psi = Omega everywhere.
%
% This serves as a direct comparison partner for SRDKF: same stochastic
% trigger, same network fusion, but no model-uncertainty robustness.
%
classdef SDKF
  properties
    Ts
    T
    % Tuning parameters
    alpha
    beta
    delta
    Z
    % Network
    G
    W
    N
    S
    % State estimate history
    X_hat
    % Model matrices
    A
    C
    Q
    R
    n
    % Stats
    RMSE
    txRate
    % 'open' or 'closed'
    triggerMode
  end

  methods
    function self = SDKF(plant, Ts, T, G, alpha, beta, delta, Z, triggerMode)
      self.Ts = Ts;
      self.T  = T;
      self.G  = G;

      self.alpha = alpha;
      self.beta  = beta;
      self.delta = delta;
      self.Z     = Z;

      self.N = numnodes(G);
      self.S = sum(G.Nodes.isSensor);

      self.W = calcMetropolisWeights(G);

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      self.X_hat  = zeros(self.n, self.N, T + 1);
      self.txRate = zeros(self.T + 1, 1);

      self.triggerMode = triggerMode;
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

        c_t = self.exchange(self.X_hat(:, :, t), Omega_upd, q_bar, Omega_bar, y);
        self.txRate(t) = sum(c_t) / self.N;

        [q_fused, Omega_fused] = self.fusion(c_t, q_upd, Omega_upd, q_bar, Omega_bar);
        [q_pred, Omega_pred]   = self.getLocalPriors(q_fused, Omega_fused);
        [q_bar, Omega_bar]     = self.updateGlobalPriors(c_t, q_upd, Omega_upd, q_bar, Omega_bar);
      end

      self.RMSE = calculateRmse(self.X_hat, X);
    end

    %% Correction step
    function [q_upd, Omega_upd] = update(self, q_pred, Omega_pred, y)
      q_upd     = zeros(self.n, self.N);
      Omega_upd = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        if self.G.Nodes(i, :).isSensor
          idx = (2 * i - 1):(2 * i);
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
    % Sensor nodes use the stochastic event trigger (Han et al. 2015);
    % non-sensor relay nodes use the DKF deterministic estimate-change
    % trigger (Battistelli et al. 2018), so information still diffuses
    % through the network backbone.
    function c_t = exchange(self, X_hat, Omega_upd, q_bar, Omega_bar, y)
      c_t = ones(self.N, 1);
      if any(isnan(Omega_bar), 'all')
        return
      end

      for i = 1:self.N
        if ~self.G.Nodes(i, :).isSensor
          % Relay node: deterministic trigger on the estimate change.
          Omega_i = Omega_upd(:, :, i);
          x_bar_i = Omega_bar(:, :, i) \ q_bar(:, i);

          e     = X_hat(:, i) - x_bar_i;
          eNorm = e' * Omega_i * e;

          lowerB = (1 / (1 + self.beta)) * Omega_i;
          upperB = (1 + self.delta) * Omega_i;

          if eNorm <= self.alpha && loewnerBetweenEig(lowerB, Omega_bar(:, :, i), upperB)
            c_t(i) = 0;
          end
          continue
        end

        idx = (2 * i - 1):(2 * i);
        y_i = y(idx);

        switch lower(self.triggerMode)
          case 'closed'
            x_bar_i = Omega_bar(:, :, i) \ q_bar(:, i);
            C_i     = self.C(idx, :);
            c_t(i)  = checkClosedLoopStochasticFusionConditions(y_i, C_i, x_bar_i, self.Z);

          case 'open'
            c_t(i) = checkOpenLoopStochasticFusionConditions(y_i, self.Z);

          otherwise
            error('Unknown triggerMode. Use ''closed'' or ''open''.');
        end
      end
    end

    %% Information fusion
    function [q_fused, Omega_fused] = fusion(self, c_t, q_upd, Omega_upd, q_bar, Omega_bar)
      q_fused     = zeros(self.n, self.N);
      Omega_fused = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        [~, nids] = inedges(self.G, i);

        for j = nids'
          w_ij = self.W(i, j);

          if (i == j) || c_t(j)
            q_fused(:, i)       = q_fused(:, i)       + w_ij * q_upd(:, j);
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + w_ij * Omega_upd(:, :, j);
          else
            q_tilde     = (1 / (1 + self.delta)) * q_bar(:, j);
            Omega_tilde = (1 / (1 + self.delta)) * Omega_bar(:, :, j);

            q_fused(:, i)       = q_fused(:, i)       + w_ij * q_tilde;
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + w_ij * Omega_tilde;
          end
        end
      end
    end

    %% Prediction step (no robust treatment)
    function [q_pred, Omega_pred] = getLocalPriors(self, q_fused, Omega_fused)
      q_pred     = zeros(self.n, self.N);
      Omega_pred = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        Omega_pred(:, :, i) = self.updateOmega(Omega_fused(:, :, i));
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

        Omega_bar_next(:, :, i) = self.updateOmega(Omega_check);
        q_bar_next(:, i)        = Omega_bar_next(:, :, i) * self.A * (Omega_check \ q_check);
      end
    end

    function newOmega = updateOmega(self, Omega)
      invQ  = self.Q \ eye(self.n);
      invQA = invQ * self.A;
      foo   = (Omega + self.A' * invQA) \ eye(self.n);
      newOmega = invQ - invQA * foo * invQA';
    end

    %% Plotting
    function plotTrajectory(self, X)
      meanX_hat = squeeze(mean(self.X_hat, 2));

      figure
      plot(meanX_hat(3, :), meanX_hat(4, :));
      hold on
      plot(X(3, :), X(4, :));
      hold off
      title(sprintf('SDKF Estimated Trajectory (%s-loop)', self.triggerMode))
      xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
      ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
      legend({sprintf('SDKF-%s', self.triggerMode), 'Actual Model'})
      grid()
    end
  end
end
