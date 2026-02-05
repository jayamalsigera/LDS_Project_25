%% State Space Model for simulating the plant systems
classdef SSModel
  properties
    A
    B
    C
    D
    n
    m
    p
  end
  methods
    function self = SSModel(A, B, C, D)
      self.A = A;
      self.B = B;
      self.C = C;
      self.D = D;

      self.n = size(A, 1);
      self.m = size(B, 2);
      self.p = size(C, 1);
    end

    function x_t = stateEq(self, x_prev)
      w_t = randn(self.m, 1);
      x_t = self.A * x_prev + self.B * w_t;
    end

    function y_t = outputEq(self, x_t)
      v_t = randn(self.p, 1);
      y_t = self.C * x_t + self.D * v_t;
    end

    function [y_t, x_t] = update(self, x_prev)
      x_t = self.stateEq(x_prev);
      y_t = self.outputEq(x_t);
    end
  end
end
