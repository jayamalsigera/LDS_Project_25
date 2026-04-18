%% Solve gamma(Omega, theta) = b by bisection
%
% gamma(Omega,theta) as defined in Section 2 of Ghion–Zorzi (2023)
%
% gamma(Omega, theta) = 1/2 * { trace[(I - theta*Omega^{-1})^{-1} - I]
%                              + log det(I - theta*Omega^{-1}) }

function theta = findTheta(Omega, b)

  % No model uncertainty
  if b <= 0
    theta = 0;
    return
  end

  n = size(Omega, 1);
  I = eye(n);

  % theta must satisfy 0 < theta < min(eig(Omega))
  theta_low = 0;
  theta_high = min(eig(Omega)) - 1e-10;

  % Bisection
  tol = 1e-9;
  max_iter = 200;

  for k = 1:max_iter
    theta = 0.5 * (theta_low + theta_high);

    M = I - theta * inv(Omega);
    gamma_val = 0.5 * (trace(inv(M) - I) + log(det(M)));

    if abs(gamma_val - b) < tol
      return
    end

    if gamma_val < b
      theta_low = theta;
    else
      theta_high = theta;
    end
  end

  theta = 0.5 * (theta_low + theta_high);

end
