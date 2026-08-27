%% Distributed Kalman Filter
%
% Implementation based on Battistelli, et. all (2018).
%
classdef DKF
  properties
    Ts
    T
    % Tunning Parameters
    alpha
    beta
    delta
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
    Qinv
    R
    n
    senBlock
    % Stats
    RMSE
    txRate
  end

  methods
    function self = DKF(plant, Ts, T, G, alpha, beta, delta)
      self.Ts = Ts;
      self.T = T;
      self.G = G;

      self.alpha = alpha;
      self.beta = beta;
      self.delta = delta;

      self.N = numnodes(G);
      self.S = sum(G.Nodes.isSensor);

      % (Battistelly, et. all, 2018) uses Metropolis Weights, but we chose to
      % adopt fusion weights from (Ghion & Zorzi, 2023) for direct comparison.
      % self.Pi = calcMetropolisWeights(G);
      self.Pi = calcFusionWeights(G);

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.Qinv = self.Q \ eye(size(self.Q));
      self.R = plant.D * plant.D';

      self.n = plant.n;

      % Measurement rows per sensor (p_i); derived from C. See RDKF.
      assert(mod(plant.p, self.S) == 0, ...
        'DKF:senBlock', 'plant.p (%d) is not divisible by sensor count (%d).', ...
        plant.p, self.S);
      self.senBlock = plant.p / self.S;

      self.X_hat = zeros(self.n, self.N, T + 1);
      self.txRate = zeros(self.T + 1, 1);
    end

    %% Estimation Method
    function self = estimate(self, x0_hat, P0, X, Y)
      % It's unrealistic to assume all nodes share same initial conditions
      % (Battistelli & Chisci, 2014), specially with "perfect knowledge", but
      % this allow us to get results in similar shape to (Ghion & Zorzi, 2023)
      q_pred = repmat(P0 \ x0_hat, 1, self.N);
      Omega_pred = repmat(P0 \ eye(self.n), 1, 1, self.N);

      self.X_hat(:, :, 1) = repmat(x0_hat, 1, self.N);

      % Initializing the "global" predictions, assuming c_t = 1 for all nodes
      % in the first iteration (i.e. the first fusion step only relies on the
      % local filtered values).
      q_bar = nan(self.n, self.N);
      Omega_bar = nan(self.n, self.n, self.N);

      for t = 2:self.T + 1
        y = Y(:, t);
        [q_upd, Omega_upd] = self.update(q_pred, Omega_pred, y);

        for i = 1:self.N
          % Return to state representation from information form
          self.X_hat(:, i, t) = pinv(Omega_upd(:, :, i)) * q_upd(:, i);
        end

        c_t = self.exchange(self.X_hat(:, :, t), Omega_upd, q_bar, Omega_bar);
        self.txRate(t) = sum(c_t) / self.N;

        [q_fused, Omega_fused] = self.fusion(c_t, q_upd, Omega_upd, q_bar, Omega_bar);
        [q_pred, Omega_pred] = self.getLocalPriors(q_fused, Omega_fused);
        [q_bar, Omega_bar] = self.updateGlobalPriors(c_t, q_upd, Omega_upd, q_bar, Omega_bar);
      end

      self.RMSE = calculateRmse(self.X_hat, X);
    end

    %% Correction/Update/Measurement step
    % Update the local information pair to obtain q_k|k and Omega_k|k of each
    % node
    function [q_upd, Omega_upd] = update(self, q_pred, Omega_pred, y)
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
          Omega_upd(:, :, i) = Omega_pred(:, :, i) + C_i' * (R_i \ C_i);
        else
          q_upd(:, i) = q_pred(:, i);
          Omega_upd(:, :, i) = Omega_pred(:, :, i);
        end
      end
    end

    %% Information Exchange
    % While we don't have to transmit anything, in this step we're calculating
    % c^i_t for all nodes.
    function c_t = exchange(self, X_hat, Omega_upd, q_bar, Omega_bar)
      c_t = ones(self.N, 1);
      if any(isnan(Omega_bar), "all")
        return % Every node transmits on the first iteration
      end

      for i = 1:self.N
        Omega_i = Omega_upd(:, :, i);

        x_bar_i = Omega_bar(:, :, i) \ q_bar(:, i);

        e = X_hat(:, i) - x_bar_i; % Discrepancy from Prediction since last transmission
        eNorm = e' * Omega_i * e; % Weighted Euclidean Norm

        lower = (1 / (1 + self.beta)) * Omega_i;
        upper = (1 + self.delta) * Omega_i;

        if eNorm <= self.alpha && loewnerBetweenEig(lower, Omega_bar(:, :, i), upper)
          c_t(i) = 0;
        end
      end
    end

    %% Information Pair Fusion
    function [q_fused, Omega_fused] = fusion(self, c_t, q_upd, Omega_upd, q_bar, Omega_bar)
      q_fused = zeros(self.n, self.N);
      Omega_fused = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        [~, nids] = inedges(self.G, i);

        for j = nids'
          pi_ij = self.Pi(i, j);

          if (i == j) || c_t(j)
            % Node i has received from j or this is a self-loop (node has
            % access to its own local info)
            q_fused(:, i) = q_fused(:, i) + pi_ij * q_upd(:, j);
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + pi_ij * Omega_upd(:, :, j);
          else
            q_tilde = (1 / (1 + self.delta)) * q_bar(:, j);
            Omega_tilde = (1 / (1 + self.delta)) * Omega_bar(:, :, j);

            q_fused(:, i) = q_fused(:, i) + pi_ij * q_tilde;
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + pi_ij * Omega_tilde;
          end
        end
      end
    end

    %% Prediction Step
    function [q_pred, Omega_pred] = getLocalPriors(self, q_fused, Omega_fused)
      q_pred = zeros(self.n, self.N);
      Omega_pred = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        q_i_F = q_fused(:, i);
        Omega_i_F = Omega_fused(:, :, i);

        Omega_pred(:, :, i) = predictOmega(Omega_i_F, self.A, self.Qinv);

        q_pred(:, i) = Omega_pred(:, :, i) * self.A * (Omega_i_F \ q_i_F);
      end
    end

    % Propagate global priors in time when there is no transmission
    function [q_bar_next, Omega_bar_next] = updateGlobalPriors(self, c_t, q_upd, Omega_upd, q_bar, Omega_bar)
      q_bar_next = zeros(self.n, self.N);
      Omega_bar_next = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        if c_t(i)
          q_i_check = q_upd(:, i);
          Omega_i_check = Omega_upd(:, :, i);
        else
          q_i_check = q_bar(:, i);
          Omega_i_check = Omega_bar(:, :, i);
        end

        Omega_bar_next(:, :, i) = predictOmega(Omega_i_check, self.A, self.Qinv);

        q_bar_next(:, i) = Omega_bar_next(:, :, i) * self.A * (Omega_i_check \ q_i_check);
      end
    end
  end
end
