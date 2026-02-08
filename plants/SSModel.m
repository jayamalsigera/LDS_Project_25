%% State Space Model for simulating the plant systems
classdef SSModel < handle
  properties
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
    function self = SSModel(A, B, C, D, T)
      self.A = A;
      self.B = B;
      self.C = C;
      self.D = D;

      self.n = size(A, 1);
      self.m = size(B, 2);
      self.p = size(C, 1);

      % Initialize state and output history
      self.T = T;
      self.X = zeros(self.n, T + 1);
      self.Y = zeros(self.p, T + 1);
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
  end
end
