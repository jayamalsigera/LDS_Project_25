%% Single Target 2D Model
%
% Based on Ghion, Zorzi (2023) but 2D like Battistelli, et al. (2018).
%
% Parameters:
% - Ts - Sampling Period
% - S - Number of Sensors
% - T - Simulation Steps
classdef SingleTarget2dModel < handle
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
    % Simulated State and Output
    T
    X
    Y
  end
  methods
    function self = SingleTarget2dModel(Ts, S, noiseStd, T)
      self.Ts = Ts;

      self.initStateEquation(Ts)
      self.initOutputEquation(S, noiseStd)

      % Initialize state and output history
      self.T = T;
      self.X = zeros(self.n, T + 1);
      self.Y = zeros(self.p, T + 1);
    end

    function initStateEquation(self, Ts)
      Ac = [zeros(2) zeros(2); eye(2) zeros(2)];
      Bc = eye(4);

      sysc = ss(Ac, Bc, eye(4), zeros(size(Bc)));
      sysd = c2d(sysc, Ts);

      self.A = sysd.A;
      self.B = sysd.B;

      self.n = size(self.A, 1);
      self.m = size(self.B, 2);
    end

    function initOutputEquation(self, S, noiseStd)
      C_i_a = [0 0 1 0; 0 0 0 0];  % Some sensors just measure p_x
      C_i_b = [0 0 0 0; 0 0 0 1];  % Some sensors just measure p_y
      self.C = [repmat(C_i_a, S/2, 1); repmat(C_i_b, S - S/2, 1)];

      D_i = noiseStd^2;

      self.D = D_i * eye(2 * S);

      self.p = size(self.C, 1);
    end

    function x_t = stateEq(self, x_prev)
      w_t = randn(self.m, 1);
      x_t = self.A * x_prev + self.B * w_t;
    end

    function y_t = outputEq(self, x_t)
      v_t = randn(self.p, 1);
      y_t = self.C * x_t + self.D * v_t;
    end

    function simulate(self, x0)
      self.X(:, 1) = x0;
      self.Y(:, 1) = self.outputEq(x0);
      for t = 2:self.T + 1

        x_prev = self.X(:, t - 1);
        self.X(:, t) = self.stateEq(x_prev);
        self.Y(:, t) = self.outputEq(self.X(:, t));
      end
    end

    function plotTrajectory(self)
      figure
      plot(self.X(3, :), self.X(4, :));
      title("Simulated Trajectory (State)")
      xlabel('$p_x$', 'Interpreter', 'latex');
      ylabel('$p_y$', 'Interpreter', 'latex');
      grid()
    end

    function plotOutputs(self)
      figure
      t = (0:self.T) * self.Ts;
      subplot(2, 1, 1)
      plot(t, self.Y(1:2, :)');
      title("Simulated Output of Sensor 1")
      subplot(2, 1, 2)
      plot(t, self.Y(39:40, :)');
      title("Simulated Output of Sensor 20")
      xlabel('Time (s)');
    end
  end
end
