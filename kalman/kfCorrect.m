%% Kalman Filter Correction Step in Information Form
% Implementation following eq. 2 and 3 in:
%   https://isas.iar.kit.edu/pdf/Fusion17_Pfaff-IDKF.pdf
function [q_new, Omega_new] = kfCorrect(q_bar, Omega_bar, C_hat, R, y)
  q_new = q_bar + C_hat' * (R \ y');
  Omega_new = Omega_bar + C_hat' * (R \ C_hat);
end
