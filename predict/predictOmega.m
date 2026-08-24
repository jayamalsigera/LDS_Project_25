%% Predict nominal information matrix
% Based on eqs 12 and 14 of (Battistelli et al., 2018)
%
% In standard Kalman, the covariance prediction step is P_pred = APA'  + Q. In
% the information form, this becomes Omega_out = (A * Omega_in^-1 * A' + Q)^-1.
% To compute it more effiently, we use the Woodbury matrix identity.
%
function Omega_out = predictOmega(Omega_in, A, Q)
  Qinv = Q \ eye(size(Q));
  Omega_out = Qinv - Qinv * A / (A' * Qinv * A + Omega_in) * A' * Qinv;
end
