%% CLSET-KF standalone unit test — reproduces Han et al. (2015) Fig. 6/7
%
% Part B step 1 (see sdkf-fix-handoff.md §7): before grafting the
% negative-information (enlarged-noise) update into the distributed
% SDKF/SRDKF, verify a *standalone* CLSET-KF against the paper's own
% target-tracking numbers.
%
% Reference: D. Han, Y. Mo, J. Wu, S. Weerakkody, B. Sinopoli, L. Shi,
% "Stochastic Event-Triggered Sensor Schedule for Remote State Estimation,"
% IEEE TAC 60(10), 2015. Theorem 2 (CLSET-KF), eqs (22)-(26); Section VI-C,
% Figs 6-7.
%
% Singer maneuvering-target model [ref 37]:
%   x = [position; speed; acceleration]
%   A = [1 T T^2; 0 1 T; 0 0 1]
%   Q = 2*alpha*sigma_m^2 * [T^5/20 T^4/8 T^3/6; T^4/8 T^3/3 T^2/2; T^3/6 T^2/2 T]
%   C = I3, R = I3,  T=1, alpha=0.01, sigma_m^2=5.
%
% CLSET-KF recursion (per step k), gamma_k in {0,1} the transmit decision:
%   x_pred = A x;               P_pred = A P A' + Q
%   z      = y - C x_pred
%   nu     = exp(-1/2 z' Z z);  gamma = (rand > nu)   % 1 = transmit
%   K      = P_pred C' / (C P_pred C' + R + (1-gamma) inv(Z))
%   x      = x_pred + gamma * K z          % mean moves only when it transmits
%   P      = P_pred - K C P_pred           % covariance shrinks even when silent
%
% Validation logic: Theorem 2 states the filter-reported P_k IS the exact
% conditional error covariance. So the mean of the filter's own P_k across
% MC runs ("theoretical") must match the sample error covariance from the
% actual estimation errors ("empirical"), and both must hit the paper's
% asymptotic P11.
%
% As a foil we also run a NAIVE filter with the *same* stochastic trigger
% but which discards the information in silence: when it does not transmit
% it performs no measurement update at all (x=x_pred, P=P_pred). This is the
% single-sensor analogue of Cause B (see handoff §5): the current SDKF/SRDKF
% reuse a discounted stale prior instead of the enlarged-noise update, so
% silence carries no information. The naive filter shares CLSET-KF's exact
% transmit decisions, so any gap between them is attributable solely to the
% negative-information update -- which is precisely what Part B adds.
%
% NOTE: this is deliberately NOT the paper's DET-KF (Wu et al. [25]). That
% scheme uses an approximate MMSE with numerical integration; reproducing its
% specific P11 numbers is out of scope for the fix. The naive foil here is a
% cleaner, more directly relevant demonstration for our distributed bug.

clear; clc;
rng(1);

%% Model
T = 1; alpha = 0.01; sigma_m2 = 5;
A = [1 T T^2; 0 1 T; 0 0 1];
C = eye(3);
R = eye(3);
Q = 2 * alpha * sigma_m2 * [T^5/20, T^4/8, T^3/6;
                            T^4/8,  T^3/3, T^2/2;
                            T^3/6,  T^2/2, T];
Lq = chol(Q + 1e-12*eye(3), 'lower');   % sampler for process noise
Lr = chol(R, 'lower');                  % sampler for measurement noise

n = 3;
Sigma0 = eye(n);        % x0 ~ N(0, Sigma0); asymptotic P is independent of this
L0 = chol(Sigma0, 'lower');

K_STEPS  = 100;
MC       = 5000;        % paper uses 10000; 5000 is enough to separate the methods
SS_WIN   = 51:K_STEPS+1; % asymptotic window (index 1 = k=0)

% (Z scale, paper empirical P11, paper theoretical P11, paper rate)
experiments = {
  0.52,  0.7991, 0.7994, 0.65;
  0.047, 4.6367, 4.6301, 0.25;
};

fprintf('CLSET-KF unit test — %d MC runs, %d steps, asymptotic window k=%d..%d\n\n', ...
        MC, K_STEPS, SS_WIN(1)-1, SS_WIN(end)-1);

for e = 1:size(experiments, 1)
  zScale   = experiments{e, 1};
  paperEmp = experiments{e, 2};
  paperThe = experiments{e, 3};
  paperRt  = experiments{e, 4};
  Z    = zScale * eye(3);
  Zinv = inv(Z);

  % Accumulators over MC runs, per time index (1 = k=0, ..., K_STEPS+1 = k=100)
  errSq_cl = zeros(K_STEPS+1, 1);   % empirical: mean actual (x1-xhat1)^2   (CLSET-KF)
  Prep_cl  = zeros(K_STEPS+1, 1);   % theoretical: mean filter P(1,1)        (CLSET-KF)
  errSq_nv = zeros(K_STEPS+1, 1);   % empirical error^2                      (naive: silence discarded)
  Prep_nv  = zeros(K_STEPS+1, 1);   % filter-reported P(1,1)                 (naive)
  txCount  = 0;                     % transmissions (shared trigger), for the rate

  for m = 1:MC
    % --- truth init and both filters' init ---
    x   = L0 * randn(n, 1);
    xc  = zeros(n, 1);  Pc = Sigma0;   % CLSET-KF (enlarged-noise silent update)
    xv  = zeros(n, 1);  Pv = Sigma0;   % naive (silence discarded)
    errSq_cl(1) = errSq_cl(1) + (x(1) - xc(1))^2;   Prep_cl(1) = Prep_cl(1) + Pc(1,1);
    errSq_nv(1) = errSq_nv(1) + (x(1) - xv(1))^2;   Prep_nv(1) = Prep_nv(1) + Pv(1,1);

    for k = 2:K_STEPS+1
      % truth + measurement
      x = A * x + Lq * randn(n, 1);
      y = C * x + Lr * randn(n, 1);

      % ---- shared stochastic transmit decision (based on CLSET-KF's innovation) ----
      xcp = A * xc;   Pcp = A * Pc * A' + Q;
      z   = y - C * xcp;
      nu  = exp(-0.5 * (z' * Z * z));
      gam = double(rand > nu);            % 1 = transmit
      txCount = txCount + gam;

      % ---- CLSET-KF: enlarged-noise update; covariance shrinks even when silent ----
      Kc  = (Pcp * C') / (C * Pcp * C' + R + (1 - gam) * Zinv);
      xc  = xcp + gam * (Kc * z);
      Pc  = Pcp - Kc * C * Pcp;
      errSq_cl(k) = errSq_cl(k) + (x(1) - xc(1))^2;
      Prep_cl(k)  = Prep_cl(k)  + Pc(1,1);

      % ---- naive: same gamma, but silence carries no information ----
      xvp = A * xv;   Pvp = A * Pv * A' + Q;
      if gam
        zv = y - C * xvp;
        Kv = (Pvp * C') / (C * Pvp * C' + R);
        xv = xvp + Kv * zv;
        Pv = Pvp - Kv * C * Pvp;
      else
        xv = xvp;      % no update, and P NOT tightened -> over-optimistic when it does shrink,
        Pv = Pvp;      % over-pessimistic in state error; the point is empirical != reported
      end
      errSq_nv(k) = errSq_nv(k) + (x(1) - xv(1))^2;
      Prep_nv(k)  = Prep_nv(k)  + Pv(1,1);
    end
  end

  errSq_cl = errSq_cl / MC;  Prep_cl = Prep_cl / MC;
  errSq_nv = errSq_nv / MC;  Prep_nv = Prep_nv / MC;

  emp_cl = mean(errSq_cl(SS_WIN));   the_cl = mean(Prep_cl(SS_WIN));
  emp_nv = mean(errSq_nv(SS_WIN));   the_nv = mean(Prep_nv(SS_WIN));
  rate   = txCount / (MC * K_STEPS);

  dev_cl = 100 * abs(emp_cl - the_cl) / the_cl;   % empirical vs theoretical (our run)
  dev_nv = 100 * abs(emp_nv - the_nv) / the_nv;
  err_cl = 100 * abs(emp_cl - paperEmp) / paperEmp; % empirical vs paper

  fprintf('=== Experiment %d:  Z = %.3f*I3  (target rate %.2f) ===\n', e, zScale, paperRt);
  fprintf('  CLSET-KF   rate=%.3f   asymptotic P11: empirical=%.4f  theoretical=%.4f  (emp-vs-theo dev=%.3f%%)\n', ...
          rate, emp_cl, the_cl, dev_cl);
  fprintf('             paper:      rate=%.2f          empirical=%.4f  theoretical=%.4f  (ours vs paper emp=%.3f%%)\n', ...
          paperRt, paperEmp, paperThe, err_cl);
  fprintf('  naive      (same trigger, silence discarded): empirical=%.4f  theoretical=%.4f  (emp-vs-theo dev=%.3f%%)\n', ...
          emp_nv, the_nv, dev_nv);
  fprintf('\n');
end

fprintf('=== CLSET-KF UNIT TEST DONE ===\n');
