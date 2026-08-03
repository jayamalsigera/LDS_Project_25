%% Centralized Robust Kalman Filter (CRKF)
%
% Centralized (single-node) version of the robust minimax filter from
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
    R
    n
    % State estimate history
    x_hat
    % Stats
    RMSE
  end

  methods
    function self = CRKF(plant, Ts, T, b)
      self.Ts = Ts;
      self.T  = T;
      self.b  = b;

      self.A = plant.A;
      self.C = plant.C;
      self.Q = plant.B * plant.B';
      self.R = plant.D * plant.D';

      self.n = plant.n;

      self.x_hat = zeros(self.n, T + 1);
    end

    %% Estimation
    function self = estimate(self, x0_hat, P0, X, Y)
      Omega = P0 \ eye(self.n);
      q     = Omega * x0_hat;

      self.x_hat(:, 1) = x0_hat;

      for t = 2:self.T + 1
        % Prediction (robust)
        [q, Omega] = self.prediction(q, Omega);

        % Correction
        y = Y(:, t);
        [q, Omega] = self.correction(q, Omega, y);

        self.x_hat(:, t) = Omega \ q;
      end

      self.RMSE = calculateRmse(reshape(self.x_hat, self.n, 1, []), X);
    end

    %% Robust prediction step — eq. (8) in Ghion & Zorzi (2023)
    function [q_pred, Psi_pred] = prediction(self, q, Omega)
      Qinv = self.Q \ eye(self.n);

      % Nominal information prediction
      Omega_pred = Qinv - Qinv * self.A / (self.A' * Qinv * self.A + Omega) * self.A' * Qinv;

      % Risk-sensitivity parameter
      theta = findTheta(Omega_pred, self.b);

      % Least-favorable information matrix and vector
      Psi_pred = Omega_pred - theta * eye(self.n);
      q_pred   = Psi_pred * self.A * (Omega \ q);
    end

    %% Correction step (same as CKF)
    function [q_upd, Omega_upd] = correction(self, q_pred, Omega_pred, y)
      Omega_upd = Omega_pred + self.C' * (self.R \ self.C);
      q_upd     = q_pred     + self.C' * (self.R \ y);
    end
  end
end
