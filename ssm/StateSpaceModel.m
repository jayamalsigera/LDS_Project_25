%% State Space Model CLass
% Represent a generic LTI SSM and compute state and output updates in steps.
classdef StateSpaceModel
  properties
    A
    B
    C
    D
  end
  methods
    function self = StateSpaceModel(A, B, C, D)
      self.A = A;
      self.B = B;
      self.C = C;
      self.D = D;
    end

    % TODO: Does it make sense to manage the state outside of the class?
    function [y, x_new] = update(self, x)
      % Generate Noise
      w = normrnd(0, 1);
      v = normrnd(0, 1);

      % Update state and output
      x_new = self.A * x + self.B * w;

      % TODO: Shouldn't `y` be computed from `x`, not `x_new`?
      y = self.C * x_new + self.D * v;
    end
  end
end