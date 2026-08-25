%% Robust prediction and propagation step for one node
%
% This function serves both eqs. (12) and (14) in (Ghion & Zorzi, 2023).
%
function [q_out, Psi_out, Omega_out, theta] = doRobustPrediction(q_in, Omega_in, A, Q, b)

  % Nominal information matrix
  Omega_out = predictOmega(Omega_in, A, Q);

  % Find theta such that gamma(Omega_out, theta) = b
  theta = findTheta(Omega_out, b);

  % Least-favorable information matrix
  Psi_out = Omega_out - theta * eye(size(Omega_out));

  % Robust predicted information vector
  q_out = Psi_out * A * (Omega_in \ q_in);
end
