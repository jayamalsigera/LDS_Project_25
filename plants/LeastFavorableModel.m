%% Least-Favorable Model (Levy & Nikoukhah 2013, Section V)
%
% Builds the worst-case state-space model in the KL ball of radius b around a
% nominal plant, and generates realizations from it. Ghion & Zorzi (2023) use
% this model to generate Monte Carlo trajectories in their simulations.
%
% The LFM is the 2n-dimensional linear Gaussian system (Eq. 94 of the paper):
%
%   xi_{t+1} = Atil_t xi_t + Btil_t eps_t
%   y_t      = Ctil_t xi_t + Dtil_t eps_t
%
% with xi_t = [x_t; e_t], e_t = x_t - x_hat_t, eps_t ~ N(0, I).
%
% Equivalently, sampling v_t = H_t e_t + L_t eps_t and driving the nominal
% plant with v_t gives the same distribution; this implementation uses that
% equivalent form for clarity.
%
% G_t (Kalman gain, Eq. 71), H_t (Eq. 87), L_t (Eq. 92) come from a
% centralized risk-sensitive Kalman sweep. The risk-sensitivity parameter
% lambda_t is the reciprocal of the theta returned by findTheta.
%
% Usage:
%   lfm = LeastFavorableModel(plant, P0, b, T);
%   lfm = lfm.simulate(x0);      % populates lfm.X, lfm.Y
%
classdef LeastFavorableModel
  properties
    Ts
    T
    b
    % Nominal matrices (with process/measurement drivers stacked into one
    % m = plant.m + plant.p dimensional noise eps_t ~ N(0, I_m)).
    A
    C
    Btil   % [plant.B, 0]     (n x m)
    Dtil   % [0, plant.D]     (p x m)
    n
    p
    m
    % Per-step LF statistics (indexed 1..T+1, corresponding to t = 0..T)
    G       % n x p x (T+1)   Kalman gain
    V       % n x n x (T+1)   LF error covariance
    lambda  % (T+1) x 1       per-step Lagrange multiplier
    H       % m x n x (T+1)   LF noise mean shift
    L       % m x m x (T+1)   LF noise covariance Cholesky factor
    % Outputs
    X
    Y
  end

  methods
    function self = LeastFavorableModel(plant, P0, b, T)
      if b <= 0
        error('LeastFavorableModel:NonPositiveTolerance', ...
          'b must be > 0; for b = 0 the LFM coincides with the nominal plant.');
      end

      self.Ts = plant.Ts;
      self.T  = T;
      self.b  = b;

      self.A = plant.A;
      self.C = plant.C;
      self.n = plant.n;
      self.p = plant.p;

      self.Btil = [plant.B, zeros(self.n, self.p)];
      self.Dtil = [zeros(self.p, plant.m), plant.D];
      self.m    = size(self.Btil, 2);

      [self.G, self.V, self.lambda, self.H, self.L] = self.precompute(P0);

      self.X = zeros(self.n, T + 1);
      self.Y = zeros(self.p, T + 1);
    end

    %% Forward risk-sensitive Kalman sweep + backward W-recursion
    function [G, V, lambda, H, L] = precompute(self, P0)
      n_ = self.n; p_ = self.p; m_ = self.m; T_ = self.T;
      A_ = self.A; C_ = self.C; Bt_ = self.Btil; Dt_ = self.Dtil;
      In = eye(n_); Im = eye(m_);

      G      = zeros(n_, p_, T_ + 1);
      V      = zeros(n_, n_, T_ + 1);
      lambda = zeros(T_ + 1, 1);
      Pnext  = zeros(n_, n_, T_ + 1);

      V(:, :, 1) = P0;

      for idx = 1:T_ + 1
        Vt = V(:, :, idx);

        % Nominal joint covariance blocks for z_t = (x_{t+1}, y_t) given Y_{t-1}
        Py  = C_ * Vt * C_' + Dt_ * Dt_';       Py  = symm(Py);
        Pxy = A_ * Vt * C_' + Bt_ * Dt_';
        Px  = A_ * Vt * A_' + Bt_ * Bt_';       Px  = symm(Px);

        Gt  = Pxy / Py;                          % Eq. 71
        Pn  = symm(Px - Pxy * (Py \ Pxy'));       % nominal one-step error cov

        G(:, :, idx)     = Gt;
        Pnext(:, :, idx) = Pn;

        % gamma(Omega, theta) = b, with Omega = Pn^{-1}, theta = 1/lambda
        Omega = Pn \ In;
        theta = findTheta(Omega, self.b);
        if theta <= 0
          error('LeastFavorableModel:InvalidTheta', ...
            'findTheta failed at t=%d (theta=%.3g); try a smaller b.', idx-1, theta);
        end
        lambda(idx) = 1 / theta;

        if idx < T_ + 1
          % Eq. 74: V_{t+1}^{-1} = P_{t+1}^{-1} - lambda^{-1} I
          invV = symm(Omega - theta * In);
          V(:, :, idx + 1) = symm(invV \ In);
        end
      end

      H = zeros(m_, n_, T_ + 1);
      L = zeros(m_, m_, T_ + 1);

      % Backward recursion on W (n x n). Boundary (Eq. 90): W_{T+1} = lambda_T I.
      W = lambda(T_ + 1) * In;

      for idx = T_ + 1:-1:1
        Gt   = G(:, :, idx);
        BtGD = Bt_ - Gt * Dt_;                   % n x m
        AtGC = A_ - Gt * C_;                     % n x n

        KvInv = symm(Im - BtGD' * (W \ BtGD));    % Eq. 86
        Kv    = symm(KvInv \ Im);

        H(:, :, idx) = Kv * BtGD' * (W \ AtGC);  % Eq. 87
        % Eq. 92: L_t L_t' = Kv. Small jitter guards against numerical slip.
        L(:, :, idx) = chol(Kv + 1e-12 * Im, 'lower');

        if idx > 1
          % Eq. 89 + Eq. 83: W_t = (Omega_t^{-1} + lambda_{t-1}^{-1} I)^{-1}
          M        = symm(W - BtGD * BtGD');
          OmegaInv = symm(AtGC' * (M \ AtGC));
          W        = symm((OmegaInv + (1 / lambda(idx - 1)) * In) \ In);
        end
      end
    end

    %% Forward-simulate the LFM to produce (X, Y)
    function self = simulate(self, x0)
      n_ = self.n; m_ = self.m; T_ = self.T;
      A_ = self.A; C_ = self.C; Bt_ = self.Btil; Dt_ = self.Dtil;

      xi = zeros(2 * n_, T_ + 1);
      xi(1:n_, 1) = x0;          % e_0 = 0

      for idx = 1:T_ + 1
        x_t = xi(1:n_, idx);
        e_t = xi(n_+1:2*n_, idx);

        Ht = self.H(:, :, idx);
        Lt = self.L(:, :, idx);
        Gt = self.G(:, :, idx);

        eps = randn(m_, 1);
        v   = Ht * e_t + Lt * eps;  % LF noise: mean H e, cov L L' = Kv

        self.X(:, idx) = x_t;
        self.Y(:, idx) = C_ * x_t + Dt_ * v;

        if idx < T_ + 1
          xi(1:n_,     idx + 1) = A_ * x_t + Bt_ * v;
          xi(n_+1:2*n_, idx + 1) = (A_ - Gt * C_) * e_t + (Bt_ - Gt * Dt_) * v;
        end
      end
    end
  end
end

function M = symm(M)
  M = (M + M') / 2;
end
