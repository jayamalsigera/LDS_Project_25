%% Distributed State Estimation Algorithm with Consensus on the Posteriors (DSEA-CP)
classdef DSEACP
  properties
    Ts
    T
    % Network Graph and Parameters
    G
    L
    % Information Pair history
    q % n x (T+1)
    Omega % n x n x (T+1)
    % State estimate history
    x_hat % n x (T+1)
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

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      self.Omega = zeros(self.n, self.n, T + 1);
      self.q = zeros(self.n, T + 1);
      self.x_hat = zeros(self.n, T + 1);
    end

    function self = estimate(self, x0_hat, P0, X, Y)
      self.Omega(:, :, 1) = inv(P0);
      self.q(:, 1) = self.Omega(:, :, 1) * x0_hat;
      self.x_hat(:, 1) = x0_hat;

      for t = 2:self.T + 1
        t_prev = t - 1; % index of previous filtered

        y = Y(:, t);

        self.x_hat(:, t) = Omega_upd \ q_upd;
        self.q(:, t) = q_upd;
        self.Omega(:, :, t) = Omega_upd;
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
