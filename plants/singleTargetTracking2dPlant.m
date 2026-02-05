%% Single Target Tracking Model
%
% Based on the simulations of the paper Battistelli, Chisci, Selvi, 2018.
%
function [plant] = singleTargetTracking2dPlant(alpha, sigma)
  Ts = 1;  % Sampling Period in seconds

  Ac_block = [0 1; 0 -alpha]
  Ac = [Ac_block zeros(2); zeros(2) Ac_block];
  Bc = [0 0; sigma 0; 0 0; 0 sigma];

  sysc = ss(Ac, Bc, eye(4), zeros(4));
  sysd = c2d(sysc, Ts);

  A = sysd.A;
  B = sysd.B;
  C_i = [1 0 0 0; 0 0 1 0];
  D_i = 10^2;  % std of 10m

  % TODO: Actually stack this for how many sensors we'll have
  C = C_i;
  % std of each component of the measurement is set to 10m
  D = D_i;

  plant = SSModel(A, B, C, D);
end
