%% Returning from information form to statespace-form state vectors
function [X, x_bar] = kfGetState(Omega, q, Omega_bar, q_bar)
  % Modification of first equation in III B in:
  % https://isas.iar.kit.edu/pdf/Fusion17_Pfaff-IDKF.pdf
  X = Omega \ q;
  x_bar = Omega_bar \ q_bar;
end
