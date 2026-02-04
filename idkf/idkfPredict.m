%% Kalman Filter Prediction Step in Information Form
% Implementation following first eq. on p. 3, col. 2 in:
%   https://isas.iar.kit.edu/pdf/Fusion17_Pfaff-IDKF.pdf
function [q_bar, Omega_bar] = idkfPredict(q_fused, Omega_fused, A_hat, Q)
  Omega_bar = inv(A_hat * (Omega_fused \ A_hat') + Q);
  q_bar = Omega_bar * A_hat * (Omega_fused \ q_fused);
end
