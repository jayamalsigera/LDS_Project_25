%% Distributed State Estimation Algorithm with Consensus on the Posteriors (DSEA-CP)
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
    % Information Pair history
    q_upd
    q_pred
    Omega_upd
    Omega_pred
    % State estimate history
    x_hat
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

      self.q_upd = zeros(self.n, self.N, T + 1);
      self.q_pred = zeros(self.n, self.N, T + 1);
      self.Omega_upd = zeros(self.n, self.n, self.N, T + 1);
      self.Omega_pred = zeros(self.n, self.n, self.N, T + 1);

      self.x_hat = zeros(self.n, self.N, T + 1);
    end

    % TODO: This is failing because it's a value class. Need to return value instead of mutating
    function update(self, t, y)
      for i = 1:self.N
        if self.G.Nodes(i, :).isSensor
          y_i = y(i:i + 1);
          C_i = self.C(i:i + 1, :);
          R_i = self.R(i:i + 1, i:i + 1); % TODO: Review if this makes sense

          self.q_upd(:, i, t) = self.q_pred(:, i, t - 1) + C_i' * (R_i \ y_i);
        else
          self.q_upd(:, i, t) = self.q_pred(:, i, t - 1);
        end

        self.Omega_upd(:, :, i, t) = self.Omega_upd(:, :, i, t - 1) + C_i' * (R_i \ C_i);
      end
    end

    %% Information Pair Fusion via Consensus Algorithm
    function [q_fused, Omega_fused] = fusion(self, t)
      q_cons = zeros(self.n, self.N, self.L);
      Omega_cons = zeros(self.n, self.n, self.N, self.L);

      q_cons(:, :, 1) = self.q_upd(:, :, t);
      Omega_cons(:, :, :, 1) = self.Omega_upd(:, :, :, t);

      for l = 2:self.L

        for i = 1:self.N
          [~, nids] = inedges(self.G, i);

          for j = nids'
            q_j = q_cons(:, j, t);
            Omega_j = Omega_cons(:, :, j, t);

            w_ij = self.W(i, j);

            q_cons(:, i, l) = q_cons(:, i, l) + w_ij * q_j;
            Omega_cons(:, :, i, l) = Omega_cons(:, :, i, l) + w_ij * Omega_j;
          end

        end

      end

      q_fused = q_cons(:, :, end);
      Omega_fused = Omega_cons(:, :, :, end);
    end

    function prediction(self, t, q_fused, Omega_fused)
      for i = 1:self.N
        % TODO: Review/Cleanup
        q_i = q_fused(:, i);
        Omega_i = Omega_fused(:, :, i);

        I = eye(self.n);

        foo = inv(Omega_i + self.A' * (self.Q \ self.A));

        bar = self.A' \ (Omega_i * self.A');

        baz = self.A' \ Omega_i;

        self.q_pred(:, i, t) = self.A' \ ((I - Omega_i * foo) * q_i);
        self.Omega_pred(:, :, i, t) = bar - baz * foo * baz';
      end

    end

    function self = estimate(self, x0_hat, X, Y)
      % Initializing from 0
      % self.q_pred(:, :, 1) =
      % self.Omega_pred(:, :, :, 1) =
      self.x_hat(:, :, 1) = repmat(x0_hat, 1, self.N);

      for t = 2:self.T + 1
        y = Y(:, t);

        disp(self.q_upd(:, :, t))
        self.update(t, y);
        disp(self.q_upd(:, :, t))
        [q_fused, Omega_fused] = self.fusion(t);
        disp(self.q_upd(:, :, t))
        self.prediction(t, q_fused, Omega_fused);

        % disp(self.Omega_upd(:, :, :, t))
        % disp(self.q_upd(:, :, t))

        % disp(self.Omega_upd)

        self.x_hat(:, t) = self.Omega_upd(:, :, :, t) \ self.q_upd(:, :, t);
      end

      self.RMSE = vecnorm(self.x_hat - X);
    end

    function plotTrajectory(self, X)
      figure
      plot(self.x_hat(3, :), self.x_hat(4, :));
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
