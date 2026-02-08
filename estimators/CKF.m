%% Centralized Kalman Filter (CKF) - Information form
classdef CKF < handle
  properties
    T
    % Information Pair history
    q % n x (T+1)
    Omega % n x n x (T+1)
    % State estimate history
    x_hat % n x (T+1)
    % Model Matrices
    plant
    Q
    R
    % Precomputed measurement info term
    CtRinvC
    n
  end

  methods
    function self = CKF(plant, T)
      self.plant = plant;
      self.T = T;

      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      self.Omega = zeros(self.n, self.n, T + 1);
      self.q = zeros(self.n, T + 1);
      self.x_hat = zeros(self.n, T + 1);

      % Precompute C'R^{-1}C (R constant)
      self.CtRinvC = plant.C' * (self.R \ plant.C);
    end

    function [q_pred, Omega_pred] = prediction(self, q_prev, Omega_prev)
      x_prev = Omega_prev \ q_prev;
      P_prev = Omega_prev \ eye(self.n);

      x_pred = self.plant.A * x_prev;
      P_pred = self.plant.A * P_prev * self.plant.A' + self.Q;

      Omega_pred = P_pred \ eye(self.n);
      q_pred = Omega_pred * x_pred;
    end

    function [q_upd, Omega_upd] = update(self, y_t, q_pred, Omega_pred)
      Omega_upd = Omega_pred + self.CtRinvC;
      q_upd = q_pred + self.plant.C' * (self.R \ y_t);
    end

    function run(self, x0_hat, P0)
      self.Omega(:, :, 1) = P0 \ eye(self.n); % inv(P0)
      self.q(:, 1) = self.Omega(:, :, 1) * x0_hat;
      self.x_hat(:, 1) = x0_hat;

      for t = 2:self.T + 1
        t_prev = t - 1; % index of previous filtered

        Omega_prev = self.Omega(:, :, t_prev);
        q_prev = self.q(:, t_prev);

        [q_pred, Omega_pred] = self.prediction(q_prev, Omega_prev);

        y = self.plant.Y(:, t);
        [q_upd, Omega_upd] = self.update(y, q_pred, Omega_pred);

        self.x_hat(:, t) = Omega_upd \ q_upd;
        self.q(:, t) = q_upd;
        self.Omega(:, :, t) = Omega_upd;
      end
    end

    function plotTrajectory(self)
      figure
      plot(self.x_hat(3, :), self.x_hat(4, :));
      hold on
      plot(self.plant.X(3, :), self.plant.X(4, :));
      hold off
      title("CKF Estimated Trajectory")
      xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
      ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
      legend({"CKF", "Actual Model"})
      grid()
    end
  end
end
