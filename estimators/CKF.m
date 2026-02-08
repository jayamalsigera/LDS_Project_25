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
    A
    C
    Q
    R
    % Precomputed measurement info term
    CtRinvC
    n
  end

  methods
    function self = CKF(plant, T)
      self.T = T;

      self.A = plant.A;
      self.C = plant.C;

      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      self.Omega = zeros(self.n, self.n, T + 1);
      self.q = zeros(self.n, T + 1);
      self.x_hat = zeros(self.n, T + 1);

      % Precompute C'R^{-1}C (R constant)
      self.CtRinvC = self.C' * (self.R \ self.C);

    end

    function [q_pred, Omega_pred] = prediction(self, q_prev, Omega_prev)
      x_prev = Omega_prev \ q_prev;
      P_prev = Omega_prev \ eye(self.n);

      x_pred = self.A * x_prev;
      P_pred = self.A * P_prev * self.A' + self.Q;

      Omega_pred = P_pred \ eye(self.n);
      q_pred = Omega_pred * x_pred;
    end

    function [x_upd, q_upd, Omega_upd] = measurement(self, y_t, q_pred, Omega_pred)
      Omega_upd = Omega_pred + self.CtRinvC;
      q_upd = q_pred + self.C' * (self.R \ y_t);

      x_upd = Omega_upd \ q_upd;
    end

    function self = run(self, x0_hat, P0, Y)
      self.Omega(:, :, 1) = P0 \ eye(self.n); % inv(P0)
      self.q(:, 1) = self.Omega(:, :, 1) * x0_hat;
      self.x_hat(:, 1) = x0_hat;

      for t = 2:self.T + 1
        t_prev = t - 1; % index of previous filtered

        Omega_prev = self.Omega(:, :, t_prev);
        q_prev = self.q(:, t_prev);

        [q_pred, Omega_pred] = self.prediction(q_prev, Omega_prev);
        [x_upd, q_upd, Omega_upd] = self.measurement(Y(:, t - 1), q_pred, Omega_pred);

        % store
        self.x_hat(:, t) = x_upd;
        self.q(:, t) = q_upd;
        self.Omega(:, :, t) = Omega_upd;
      end
    end
  end
end
