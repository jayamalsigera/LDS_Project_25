%% Solve gamma(Omega, theta) = b by bisection
%
% gamma(Omega,theta) as defined in Section 2 of Ghion–Zorzi (2023)
%
% gamma(Omega, theta) = 1/2 * { trace[(I - theta*Omega^{-1})^{-1} - I]
%                              + log det(I - theta*Omega^{-1}) }
%
function theta = findTheta(Omega, b)
  if b <= 0  % No model uncertainty, degenerate to standard Kalman equations
    theta = 0;
    return
  end

  n = size(Omega, 1);
  I = eye(n);
  invOmega = Omega \ I;

  gamma = @(theta) 0.5 * (trace(inv(I - theta * invOmega) - I) ...
                          + log(det(I - theta * invOmega)));

  tol = 1e-9;
  maxIter = 200;

  % theta must satisfy 0 < theta < min(eig(Omega))
  thetaLow = 0;
  thetaHigh = min(eig(Omega)) - tol;

  theta = bisect(gamma, b, thetaLow, thetaHigh, tol, maxIter);
end

%% Bisection algorithm for computing theta_t
%
% Returns when either |a - b| <= tol or `maxIter` iterations have elapsed.
%
% Computed according to (Zenere & Zorzi, 2018)
function x = bisect(fun, b, thetaLow, thetaHigh, tol, maxIter)
  for k = 1:maxIter
    if abs(thetaLow - thetaHigh) <= tol
      return
    end

    x = 0.5 * (thetaLow + thetaHigh);
    fx = fun(x);

    if fx < b
      thetaLow = x;
    else
      thetaHigh = x;
    end
  end

  warning('bisect:MaxIterReached', ...
          'Maximum iterations (%d) reached before convergence.', maxIter);
end
