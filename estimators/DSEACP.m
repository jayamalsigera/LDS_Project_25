%% Distributed State Estimation Algorithm with Consensus on the Posteriors (DSEA-CP)
%
% Implementation based on Battistelli & Chisci (2014).
%
classdef DSEACP
  properties
    Ts
    T
    % Network Graph and Parameters
    G
    W
    L
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
  end

  methods
    function self = DSEACP(plant, Ts, T, G, L)
      self.Ts = Ts;
      self.T = T;
      self.G = G;
      self.L = L;

      self.N = numnodes(G);
      self.S = sum(G.Nodes.isSensor);

      self.W = calcMetropolisWeights(G);

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      self.X_hat = zeros(self.n, self.N, T + 1);
    end

    %% Estimation Method
    function self = estimate(self, x0_hat, P0, X, Y)
      % It's unrealistic to assume all nodes share same initial conditions
      % (Battistelli & Chisci, 2014), specially with "perfect knowledge", but
      % this allow us to get results in similar shape to (Ghion & Zorzi, 2023)
      q_pred = repmat(P0 \ x0_hat, 1, self.N);
      Omega_pred = repmat(P0 \ eye(self.n), 1, 1, self.N);

      self.X_hat(:, :, 1) = repmat(x0_hat, 1, self.N);

      for t = 2:self.T + 1
        y = Y(:, t);

        [q_upd, Omega_upd] = self.update(q_pred, Omega_pred, y);
        [q_fused, Omega_fused] = self.fusion(q_upd, Omega_upd);
        [q_pred, Omega_pred] = self.prediction(q_fused, Omega_fused);

        % Pack `q_pred` vectors page-wise to compute state estimation for
        % each node without a for loop.
        % B = reshape(q_pred, self.n, 1, self.N);
        % % Use `pagelsqminnorm` which is similar to using the Moore-Penrose
        % % pseudo inversion but page-wise.
        % x_nodes = pagelsqminnorm(Omega_pred, B);
        % self.X_hat(:, :, t) = reshape(x_nodes, self.n, self.N);

        for i = 1:self.N
          self.X_hat(:, i, t) = pinv(Omega_pred(:, :, i)) * q_pred(:, i);
        end
      end

      self.RMSE = self.calculateRSME(self.X_hat, X);
    end

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

    %% Information Pair Fusion via Consensus Algorithm
    function [q_fused, Omega_fused] = fusion(self, q_upd, Omega_upd)
      q_cons = zeros(self.n, self.N, self.L);
      Omega_cons = zeros(self.n, self.n, self.N, self.L);

      q_cons(:, :, 1) = q_upd;
      Omega_cons(:, :, :, 1) = Omega_upd;

      for l = 2:self.L

        for i = 1:self.N
          [~, nids] = inedges(self.G, i);

          % TODO: Check if this is covering the self loop
          for j = nids'
            q_j = q_cons(:, j, l - 1);
            Omega_j = Omega_cons(:, :, j, l - 1);

            w_ij = self.W(i, j);

            q_cons(:, i, l) = q_cons(:, i, l) + w_ij * q_j;
            Omega_cons(:, :, i, l) = Omega_cons(:, :, i, l) + w_ij * Omega_j;
          end

        end

      end

      q_fused = q_cons(:, :, end);
      Omega_fused = Omega_cons(:, :, :, end);
    end

    function [q_pred, Omega_pred] = prediction(self, q_fused, Omega_fused)
      q_pred = zeros(self.n, self.N);
      Omega_pred = zeros(self.n, self.n, self.N);

      for i = 1:self.N
        q_i = q_fused(:, i);
        Omega_i = Omega_fused(:, :, i);

        % TODO: Review/Cleanup
        I = eye(self.n);
        foo = inv(Omega_i + self.A' * (self.Q \ self.A));
        bar = self.A' \ (Omega_i / self.A);
        baz = self.A' \ Omega_i;

        q_pred(:, i) = self.A' \ ((I - Omega_i * foo) * q_i);
        Omega_pred(:, :, i) = bar - baz * foo * baz';
      end

    end

    function [rmse] = calculateRSME(self, X_hat, X)
      rmse = zeros(self.T+1, 1);
      for t = 1:self.T
        err = X_hat(:, :, t) - X(:, t); % implicit expansion over N
        rmse(t) = sqrt(mean(sum(err .^ 2, 1)));
      end
    end

    function plotTrajectory(self, X)
      % TODO: Would be cool if we could plot P(t) somehow
      % TODO: Restrict axis to ranges of X

      meanX_hat = mean(self.X_hat, 2);

      figure
      plot(meanX_hat(3, :), meanX_hat(4, :));
      hold on
      plot(X(3, :), X(4, :));
      hold off
      title(sprintf("DSEA-CP (L=%d) Estimated Trajectory", self.L))
      xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
      ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
      legend({"DSEA-CP", "Actual Model"})
      grid()
    end
  end
end
