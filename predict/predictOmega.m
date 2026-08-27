%% Predict nominal information matrix
% Based on eqs 12 and 14 of (Battistelli et al., 2018)
%
% In standard Kalman, the covariance prediction step is P_pred = APA'  + Q. In
% the information form, this becomes Omega_out = (A * Omega_in^-1 * A' + Q)^-1.
% To compute it more efficiently, we use the Woodbury matrix identity.
%
% Takes Qinv rather than Q: the process noise is constant, so callers precompute
% inv(Q) once in their constructor instead of once per node per time step.
%
function Omega_out = predictOmega(Omega_in, A, Qinv)
  Omega_out = Qinv - Qinv * A / (A' * Qinv * A + Omega_in) * A' * Qinv;
end
