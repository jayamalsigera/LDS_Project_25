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
    function self = DSEACP(plant, G, L, Ts, T)
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

    function [q_pred, Omega_pred] = prediction(self, q_prev, Omega_prev)
      x_prev = Omega_prev \ q_prev;
      P_prev = Omega_prev \ eye(self.n);

      x_pred = self.A * x_prev;
      P_pred = self.A * P_prev * self.A' + self.Q;

      Omega_pred = P_pred \ eye(self.n);
      q_pred = Omega_pred * x_pred;
    end

    function [q_upd, Omega_upd] = update(self, y_t, q_pred, Omega_pred)
      Omega_upd = Omega_pred + self.CtRinvC;
      q_upd = q_pred + self.C' * (self.R \ y_t);
    end

    function self = estimate(self, x0_hat, P0, X, Y)
      self.Omega(:, :, 1) = P0 \ eye(self.n); % inv(P0)
      self.q(:, 1) = self.Omega(:, :, 1) * x0_hat;
      self.x_hat(:, 1) = x0_hat;

      for t = 2:self.T + 1
        t_prev = t - 1; % index of previous filtered

        Omega_prev = self.Omega(:, :, t_prev);
        q_prev = self.q(:, t_prev);

        [q_pred, Omega_pred] = self.prediction(q_prev, Omega_prev);

        y = Y(:, t);
        [q_upd, Omega_upd] = self.update(y, q_pred, Omega_pred);

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
      title("DSEA-CP Estimated Trajectory")
      xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
      ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
      legend({"DSEA-CP", "Actual Model"})
      grid()
    end
  end
end
