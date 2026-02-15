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
    W
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

      self.W = calcMetropolisWeights(G);

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      self.X_hat = zeros(self.n, self.N, T + 1);
      self.txRate = zeros(self.N, 1);
    end

    %% Estimation Method
    function self = estimate(self, x0_hat, P0, X, Y)
      % It's unrealistic to assume all nodes share same initial conditions
      % (Battistelli & Chisci, 2014), specially with "perfect knowledge", but

      q_pred = repmat(P0 \ x0_hat, 1, self.N);
      Omega_pred = repmat(P0 \ eye(self.n), 1, 1, self.N);

      self.X_hat(:, :, 1) = repmat(x0_hat, 1, self.N);

      % Initializing the "global" predictions, assuming c_t = 1 for all nodes in the first
      % iteration (i.e. the first fusion step only relies on the local filtered values).
      q_bar = nan(self.n, self.N);
      Omega_bar = nan(self.n, self.n, self.N);

      for t = 2:self.T + 1
        y = Y(:, t);
        [q_upd, Omega_upd] = self.update(q_pred, Omega_pred, y);

        for i = 1:self.N
          self.X_hat(:, i, t) = pinv(Omega_pred(:, :, i)) * q_pred(:, i);
        end

        c_t = self.exchange(self.X_hat(:, :, t), Omega_upd, q_bar, Omega_bar);
        self.txRate(t) = sum(c_t) / self.N;

        [q_fused, Omega_fused] = self.fusion(c_t, q_upd, Omega_upd, q_bar, Omega_bar);
        [q_pred, Omega_pred] = self.getLocalPriors(q_fused, Omega_fused);
        [q_bar, Omega_bar] = self.updateGlobalPriors(c_t, q_upd, Omega_upd, q_bar, Omega_bar);
      end

      self.RMSE = self.calculateRSME(self.X_hat, X);
    end

    %% Correction/Update/Measurement step
    % Update the local information pair to obtain q_k|k and Omega_k|k of each
    % node
    function [q_upd, Omega_upd] = update(self, q_pred, Omega_pred, y)
      q_upd = zeros(self.n, self.N);
      Omega_upd = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        if self.G.Nodes(i, :).isSensor
          idx = (2 * i - 1):(2 * i);

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
        % disp("OMEGA_BAR ALL NAN")
        return % Every node transmits on the first iteration
      end

      for i = 1:self.N
        Omega_i = Omega_upd(:, :, i);

        x_bar_i = Omega_bar(:, :, i) \ q_bar(:, i);

        e = X_hat(:, i) - x_bar_i; % Discrepancy from Prediction since last transmission
        eNorm = e' * Omega_i * e; % Weighted Euclidean Norm

        lower = (1 / (1 + self.beta)) * Omega_i;
        upper = (1 + self.delta) * Omega_i;

        if eNorm <= self.alpha && isPSD(Omega_bar(:, :, i) - lower) && isPSD(upper - Omega_bar(:, :, i))
          c_t(1) = 0;
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
          w_ij = self.W(i, j);

          if (i == j) || c_t(j)
            % disp(["a" i j double(i == j) c_t(j)])
            % Node i has received from j or this is a self-loop (node has access to its own local info)
            q_fused(:, i) = q_fused(:, i) + w_ij * q_upd(:, j);
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + w_ij * Omega_upd(:, :, j);
          else
            % disp(["b" i j (i == j) c_t(j)])
            q_tilda = (1 / (1 + self.delta)) * q_bar(:, j);
            Omega_tilda = (1 / (1 + self.delta)) * Omega_bar(:, :, j);

            q_fused(:, i) = q_fused(:, i) + w_ij * q_tilda;
            Omega_fused(:, :, i) = Omega_fused(:, :, i) + w_ij * Omega_tilda;
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

        Omega_pred(:, :, i) = self.updateOmega(Omega_i_F);

        q_pred(:, i) = Omega_pred(:, :, i) * self.A * (Omega_i_F \ q_i_F);
      end
    end

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

        Omega_bar_next(:, :, i) = self.updateOmega(Omega_i_check);

        q_bar_next(:, i) = Omega_bar_next(:, :, i) * self.A * (Omega_i_check \ q_i_check);
      end
    end

    function newOmega = updateOmega(self, Omega)
      invQ = self.Q \ eye(self.n);
      invQA = invQ * self.A;
      % TODO: Review variable name
      foo = (Omega + (self.A' * invQA)) \ eye(self.n);

      newOmega = invQ - invQA * foo * invQA'; % Assuming Q = Q'
    end

    %% RMSE Calculation
    % TODO: Maybe we can move this to `utils`?
    function [rmse] = calculateRSME(self, X_hat, X)
      rmse = zeros(self.T + 1, 1);
      for t = 1:self.T
        err = X_hat(:, :, t) - X(:, t); % implicit expansion over N
        rmse(t) = sqrt(mean(sum(err .^ 2, 1)));
      end
    end

    %% Plotting
    function plotTrajectory(self, X)
      % TODO: Would be cool if we could plot P(t) somehow
      % TODO: Restrict axis to ranges of X

      meanX_hat = mean(self.X_hat, 2);

      figure
      plot(meanX_hat(3, :), meanX_hat(4, :));
      hold on
      plot(X(3, :), X(4, :));
      hold off
      title("DKF Estimated Trajectory")
      xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
      ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
      legend({"DKF", "Actual Model"})
      grid()
    end
  end
end
