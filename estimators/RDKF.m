%% Robust Distributed Kalman Filter
%
% Implementation based on Ghion, Zorzi (2023).
%
classdef RDKF
  properties
    Ts
    T
    % Tunning Parameters
    alpha
    beta
    delta
    b
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
    function self = RDKF(plant, Ts, T, G, alpha, beta, delta, b)
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

      % Measurement rows per sensor (p_i). Homogeneous across sensors (3 in this
      % model), derived from the stacked C rather than hardcoded; sensors are the
      % first S nodes, so node i owns block i.
      assert(mod(plant.p, self.S) == 0, ...
        'RDKF:senBlock', 'plant.p (%d) is not divisible by sensor count (%d).', ...
        plant.p, self.S);
      self.senBlock = plant.p / self.S;

      % b may be a scalar (uniform tolerance, Algorithm 1) or an N-vector
      % (per-node local tolerances b^i, Algorithm 2 / RDKFLOC). Store it as an
      % N-vector so the prediction steps can index self.b(i) uniformly; a
      % scalar broadcasts to a constant vector, leaving scalar callers
      % unchanged (and b = 0 still recovers DKF per node).
      if isscalar(b), b = repmat(b, self.N, 1); end
      self.b = b;

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

        c_t = self.exchange(self.X_hat(:, :, t), Omega_upd, q_bar, Psi_bar);
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

    %% Information Exchange
    % Evaluate whether to transmit measurement
    function c_t = exchange(self, X_hat, Omega_upd, q_bar, Psi_bar)
      c_t = ones(self.N, 1);
      if any(isnan(Psi_bar), "all")
        return % Every node transmits on the first iteration
      end

      for i = 1:self.N
        Omega_i = Omega_upd(:, :, i);

        x_bar_i = Psi_bar(:, :, i) \ q_bar(:, i);

        e = X_hat(:, i) - x_bar_i; % Discrepancy from Prediction since last transmission
        eNorm = e' * Omega_i * e; % Weighted Euclidean Norm

        lower = (1 / (1 + self.beta)) * Omega_i;
        upper = (1 + self.delta) * Omega_i;

        if eNorm <= self.alpha && loewnerBetweenEig(lower, Psi_bar(:, :, i), upper) %Deterministic trigger condition
          c_t(i) = 0;
        end
      end
    end

    %% Information Pair Fusion
    function [q_fused, Omega_fused] = fusion(self, c_t, q_upd, Omega_upd, q_bar, Psi_bar)
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
            Omega_tilde = (1 / (1 + self.delta)) * Psi_bar(:, :, j);

            q_fused(:, i) = q_fused(:, i) + pi_ij * q_tilde;
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + pi_ij * Omega_tilde;
          end
        end
      end
    end

    %% Prediction Step
    function [q_pred, Psi, theta_vec] = getLocalPriors(self, q_fused, Omega_fused)
      q_pred = zeros(self.n, self.N);
      Omega_pred = zeros(self.n, self.n, self.N);
      Psi = zeros(self.n, self.n, self.N);
      theta_vec = zeros(self.N, 1);

      for i = 1:self.N
        q_i_F = q_fused(:, i);
        Omega_i_F = Omega_fused(:, :, i);

        Omega_pred(:, :, i) = self.updateOmega(Omega_i_F);

        % Risk Sensitivity Parameter
        theta = findTheta(Omega_pred(:, :, i), self.b(i));
        theta_vec(i) = theta;

        % Robust information Pair
        Psi(:, :, i) = Omega_pred(:, :, i) - theta * eye(self.n);
        q_pred(:, i) = Psi(:, :, i) * self.A * (Omega_i_F \ q_i_F);
      end
    end

    % Update the global priors used in the fusion rule when there is no transmission
    function [q_bar_next, Psi_bar_next, theta_bar_vec] = updateGlobalPriors(self, c_t, q_upd, Omega_upd, q_bar, Psi_bar)
      q_bar_next = zeros(self.n, self.N);
      Omega_bar_next = zeros(self.n, self.n, self.N);
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

        Omega_bar_next(:, :, i) = self.updateOmega(Omega_i_check);

        % Risk Sensitivity Parameter
        theta_bar = findTheta(Omega_bar_next(:, :, i), self.b(i));
        theta_bar_vec(i) = theta_bar;

        % Robust information pair
        Psi_bar_next(:, :, i) = Omega_bar_next(:, :, i) - theta_bar * eye(self.n);
        q_bar_next(:, i) = Psi_bar_next(:, :, i) * self.A * (Omega_i_check \ q_i_check);
      end
    end

    function newOmega = updateOmega(self, Omega)
      invQ = self.Q \ eye(self.n);
      invQA = invQ * self.A;
      % TODO: Review variable name
      foo = (Omega + (self.A' * invQA)) \ eye(self.n);

      newOmega = invQ - invQA * foo * invQA'; % Assuming Q = Q'
    end
  end
end
