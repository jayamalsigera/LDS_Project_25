%% Robust prediction step for one node
% Based on eq. (12) in Ghion–Zorzi (2023)

function [q_pred, Psi_pred, Omega_pred, theta] = predictRobustFusion(q_fused, Omega_fused, A, Q, b)

  % Nominal information prediction
  Omega_pred = inv(Q) - inv(Q) * A / (A' * inv(Q) * A + Omega_fused) * A' * inv(Q);

  % Find theta such that gamma(Omega_pred, theta) = b
  theta = solveTheta(Omega_pred, b);

  % Least-favorable information matrix
  Psi_pred = Omega_pred - theta * eye(size(Omega_pred));

  % Robust predicted information vector
  q_pred = Psi_pred * A * (Omega_fused \ q_fused);

end