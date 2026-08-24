%% Predict nominal information matrix
% Based on eqs 12 and 14 of (Battistelli et al., 2018)

function Omega_out = predictOmega(Omega_in, A, Q)

  Qinv = Q \ eye(size(Q));
  Omega_out = Qinv - Qinv * A / (A' * Qinv * A + Omega_in) * A' * Qinv;

end
