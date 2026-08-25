%% Centralized Kalman Filter (CKF) - Information form
classdef CKF
  properties
    Ts
    T
    % Model Matrices
    A
    C
    Q
    R
    n
    % Precomputed measurement info term
    CtRinvC
    % State estimate history
    x_hat % n x (T+1)
    % Stats
    RMSE
  end

  methods
    function self = CKF(plant, Ts, T)
      self.Ts = Ts;
      self.T = T;

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      self.x_hat = zeros(self.n, T + 1);

      % Precompute C'R^{-1}C (R constant)
      self.CtRinvC = self.C' * (self.R \ self.C);
    end

    %% Estimation Method
    function self = estimate(self, x0_hat, P0, X, Y)
        %Convert initial parameters to information form
      Omega_pred = P0 \ eye(self.n); % inv(P0)
      q_pred     = Omega_pred * x0_hat;

      for t = 1:self.T + 1
        [q_upd, Omega_upd] = self.update(q_pred, Omega_pred, Y(:, t));
        self.x_hat(:, t) = Omega_upd \ q_upd;

        [q_pred, Omega_pred] = self.prediction(q_upd, Omega_upd);
      end

      self.RMSE = calculateRmse(reshape(self.x_hat, self.n, 1, []), X);
    end

    %% Correction step
    function [q_upd, Omega_upd] = update(self, q_pred, Omega_pred, y)
      Omega_upd = Omega_pred + self.CtRinvC;
      q_upd     = q_pred     + self.C' * (self.R \ y);
    end

    %% Prediction Step
    % TODO: Maybe there's a way to do this using Woodbury's identity?
    function [q_pred, Omega_pred] = prediction(self, q_upd, Omega_upd)
      x_prev = Omega_upd \ q_upd;
      P_prev = Omega_upd \ eye(self.n);

      x_pred = self.A * x_prev;
      P_pred = self.A * P_prev * self.A' + self.Q;

      Omega_pred = P_pred \ eye(self.n);
      q_pred = Omega_pred * x_pred;
    end
  end
end
