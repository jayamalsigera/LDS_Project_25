%% Single Target 3D Model
%
% 3D target-tracking model from Ghion, Zorzi (2023), Section 6 (Eq. 29): a
% constant-velocity (triple-integrator) target observed by a network of sensors,
% each of which measures 2 of the 3 position coordinates.
%
% State ordering keeps velocities first, then positions (mirroring the 2D
% model): x = [vx vy vz px py pz]', n = 6.
%
% Parameters:
% - Ts - Sampling Period
% - S  - Number of Sensors
% - k  - Measurement-noise scale (knob; R^i = sqrt(k) * P_i R0 P_i')
% - T  - Simulation Steps
classdef SingleTarget3dModel
  properties
    Ts
    A
    B
    C
    D
    % Model Dimensions
    n
    m
    p
    % Per-sensor measurement-noise data
    R        % full block-diagonal measurement covariance (p x p)
    Perm     % 3 x 3 x S random permutations used to scramble R0 (inspection only)
    % Simulated State and Output
    T
    X
    Y
  end
  methods
    function self = SingleTarget3dModel(Ts, S, k, T)
      if nargin < 3 || isempty(k)
        k = 1;  % Default measurement-noise scale
      end
      self.Ts = Ts;

      % Constant-velocity (triple-integrator) dynamics. The paper uses the
      % explicit Euler form A = I + Ts*Phi (not exact c2d):
      %   d/dt [vx vy vz] = 0
      %   d/dt [px py pz] = [vx vy vz]
      Phi = [zeros(3, 3), zeros(3, 3);
             eye(3),      zeros(3, 3)];
      self.A = eye(6) + Ts * Phi;

      % Process noise: w_t ~ N(0, I_6), B = sqrt(0.001) I_6 => Q = B B' = 0.001 I.
      self.B = sqrt(0.001) * eye(6);

      self.n = size(self.A, 1);
      self.m = size(self.B, 2);

      % Measurement model: each sensor measures 2 of the 3 position coordinates,
      % laid out as a 3-row block with one identically-zero row (p_i = 3).
      %
      % The paper says each sensor measures "either two horizontal dimensions or
      % a combination of one horizontal dimension and the vertical dimension".
      % With x, y horizontal and z vertical, the second category has TWO
      % realizations, so the design has THREE types -- all three coordinate
      % pairs. The two matrices printed in the paper (diag(1,0,1), diag(0,1,1))
      % are one exemplar from each verbal category, not the complete set.
      %   horizontal-horizontal:  measures px, py
      %   horizontal-vertical #1: measures py, pz
      %   horizontal-vertical #2: measures px, pz
      % Sensors are split as evenly as possible over the three types, in
      % contiguous blocks (analogous to the 2D C_i_a / C_i_b split). This yields
      % balanced per-coordinate coverage; dropping one type would leave one
      % coordinate observed by twice as many sensors as the other two.
      C_types = {[zeros(3, 3), diag([1 1 0])], ...   % px, py  (row 3 zero)
                 [zeros(3, 3), diag([0 1 1])], ...   % py, pz  (row 1 zero)
                 [zeros(3, 3), diag([1 0 1])]};      % px, pz  (row 2 zero)
      nPerType = floor(S / 3) * ones(1, 3);
      nPerType(1:(S - sum(nPerType))) = nPerType(1:(S - sum(nPerType))) + 1;
      Cblocks = arrayfun(@(t) repmat(C_types{t}, nPerType(t), 1), 1:3, ...
                         'UniformOutput', false);
      self.C = vertcat(Cblocks{:});
      self.p = size(self.C, 1);

      % Measurement noise: R^i = sqrt(k) * P_i R0 P_i', with R0 a fixed
      % ill-conditioned base and P_i a random permutation drawn once per sensor.
      % The draws are recorded in self.Perm for inspection only -- the
      % constructor always redraws, so reproducing a given R^i requires seeding
      % the RNG (`rng`) before construction, not passing Perm back in.
      % D^i = chol(R^i)' so the additive form y = C x + D v, v ~ N(0,I),
      % reproduces R^i; the assembled D/R are block-diagonal.
      R0 = 0.5 * diag([1 4 7]);
      self.Perm = zeros(3, 3, S);
      Dblocks = cell(S, 1);
      Rblocks = cell(S, 1);
      for i = 1:S
        Pi = eye(3);
        Pi = Pi(randperm(3), :);
        self.Perm(:, :, i) = Pi;

        Ri = sqrt(k) * (Pi * R0 * Pi');
        Rblocks{i} = Ri;
        Dblocks{i} = chol(Ri, 'lower');  % D_i D_i' = R_i
      end
      self.D = blkdiag(Dblocks{:});
      self.R = blkdiag(Rblocks{:});

      % Initialize state and output history
      self.T = T;
      self.X = zeros(self.n, T + 1);
      self.Y = zeros(self.p, T + 1);
    end

    %% Update state
    function x_t = stateEq(self, x_prev)
      w_t = randn(self.m, 1);
      x_t = self.A * x_prev + self.B * w_t;
    end

    %% Generate Output
    function y_t = outputEq(self, x_t)
      v_t = randn(self.p, 1);
      y_t = self.C * x_t + self.D * v_t;
    end

    %% Iterate for all times
    function self = simulate(self, x0)
      self.X(:, 1) = x0;
      self.Y(:, 1) = self.outputEq(x0);
      for t = 2:self.T + 1
        x_prev = self.X(:, t - 1);
        self.X(:, t) = self.stateEq(x_prev);
        self.Y(:, t) = self.outputEq(self.X(:, t));
      end
    end

    %% Plot the trajectory of the states
    function plotTrajectory(self)
      figure
      plot3(self.X(4, :), self.X(5, :), self.X(6, :));
      title("Simulated Trajectory (State)")
      xlabel('$p_x$', 'Interpreter', 'latex');
      ylabel('$p_y$', 'Interpreter', 'latex');
      zlabel('$p_z$', 'Interpreter', 'latex');
      grid()
    end

    function plotOutputs(self)
      figure
      t = (0:self.T) * self.Ts;
      subplot(2, 1, 1)
      plot(t, self.Y(1:3, :)');
      title("Simulated Output of Sensor 1")
      subplot(2, 1, 2)
      plot(t, self.Y(end-2:end, :)');
      title("Simulated Output of Last Sensor")
      xlabel('Time (s)');
    end
  end
end
