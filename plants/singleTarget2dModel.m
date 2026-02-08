%% Single Target 2D Model
%
% Based on Ghion, Zorzi (2023) but 2D like Battistelli, et al. (2018).
%
% Parameters:
% - Ts - Sampling Period
% - S - Number of Sensors
% - T - Simulation Steps
%
% Returns:
% - plant - The state space model representing the full plant
%
function [plant] = singleTarget2dModel(Ts, S, T)
  Ac = [zeros(2) zeros(2); eye(2) zeros(2)];
  Bc = eye(4);

  sysc = ss(Ac, Bc, eye(4), zeros(size(Bc)));
  sysd = c2d(sysc, Ts);

  C_i_a = [0 0 1 0; 0 0 0 0];  % Some sensors just measure p_x
  C_i_b = [0 0 0 0; 0 0 0 1];  % Some sensors just measure p_y

  D_i = 10^2;  % std of 10m

  C = [repmat(C_i_a, S/2, 1); repmat(C_i_b, S - S/2, 1)];
  D = D_i * eye(2 * S);

  plant = SSModel(sysd.A, sysd.B, C, D, T);
end
