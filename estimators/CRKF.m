%% Centralized Robust Kalman Filter (CRKF)
%
% Centralized (single estimator) version of the robust minimax filter from
% Ghion & Zorzi (2023). All sensor outputs are processed together, so this
% serves as a sanity check and performance upper bound for the distributed
% RDKF / SRDKF: a well-tuned RDKF should approach but never exceed CRKF.
%
% The correction step is identical to the standard CKF. The prediction step
% is robust: it finds the risk-sensitivity parameter theta_t > 0 such that
% gamma(Omega_pred, theta_t) = b (KL tolerance), then applies
%   Psi_pred = Omega_pred - theta_t * I
%   q_pred   = Psi_pred * A * Omega^{-1} * q
%
classdef CRKF
  properties
    Ts
    T
    b
    % Model matrices
    A
    C
    Q
    Qinv
    R
    n
    % Precomputed measurement info term
    CtRinvC
    % State estimate history
    x_hat
    % Stats
    RMSE
    theta_hist
  end

  methods
    function self = CRKF(plant, Ts, T, b)
      self.Ts = Ts;
      self.T  = T;
      self.b  = b;

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.Qinv = self.Q \ eye(size(self.Q));
      self.R = plant.D * plant.D';

      self.n = plant.n;

      self.x_hat = zeros(self.n, T + 1);
      self.theta_hist = zeros(1, T + 1);

      % Precompute C'R^{-1}C (R constant)
      self.CtRinvC = self.C' * (self.R \ self.C);
    end

    %% Estimation
    function self = estimate(self, x0_hat, P0, X, Y)
      Omega_pred = P0 \ eye(self.n);
      q_pred     = Omega_pred * x0_hat;

      for t = 1:self.T + 1
        [q_upd, Omega_upd] = self.update(q_pred, Omega_pred, Y(:, t));
        self.x_hat(:, t) = Omega_upd \ q_upd;

        [q_pred, Omega_pred, theta_t] = self.prediction(q_upd, Omega_upd);
        self.theta_hist(t) = theta_t;
      end

      self.RMSE = calculateRmse(reshape(self.x_hat, self.n, 1, []), X);
    end

    %% Correction step
    function [q_upd, Omega_upd] = update(self, q_pred, Omega_pred, y)
      Omega_upd = Omega_pred + self.CtRinvC;
      q_upd     = q_pred     + self.C' * (self.R \ y);
    end

    %% Robust Prediction according to (Ghion & Zorzi, 2023)
    function [q_pred, Psi_pred, theta] = prediction(self, q_upd, Omega_upd)
      [q_pred, Psi_pred, ~, theta] = doRobustPrediction(q_upd, Omega_upd, self.A, self.Qinv, self.b);
    end
  end
end
