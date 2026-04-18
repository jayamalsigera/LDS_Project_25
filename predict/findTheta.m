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
  invOmega = Omega \ I;

  gamma = @(theta) 0.5 * (trace(inv(I - theta * invOmega) - I) ...
                          + log(det(I - theta * invOmega)));

  % theta must satisfy 0 < theta < min(eig(Omega))
  theta_low = 0;
  theta_high = min(eig(Omega)) - 1e-10;

  theta = bisect(@(th) gamma(th) - b, theta_low, theta_high, 1e-9, 200);

end
