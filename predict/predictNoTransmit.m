%% Robust propagation step for no-transmission pair
% Based on eq. (14) in Ghion–Zorzi (2023)

function [q_bar, Psi_bar, Omega_bar, theta_bar] = predictNoTransmit(q_check, Omega_check, A, Q, b)

  Qinv = Q \ eye(size(Q));

  % Nominal propagated information matrix
  Omega_bar = Qinv - Qinv * A / (A' * Qinv * A + Omega_check) * A' * Qinv;

  % Find theta_bar such that gamma(Omega_bar, theta_bar) = b
  theta_bar = findTheta(Omega_bar, b);

  % Least-favorable propagated information matrix
  Psi_bar = Omega_bar - theta_bar * eye(size(Omega_bar));

  % Propagated information vector
  q_bar = Psi_bar * A * (Omega_check \ q_check);

end
