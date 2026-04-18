function c = checkOpenLoopStochasticFusionConditions(y, Z)

  % Implement stochastic triggering using measurement directly
  % nu = exp( -1/2 * y' * Z * y )
  nu = exp(-1/2 * y' * Z * y);

  % Generate random number
  zeta = rand;

  % Transmission decision
  if zeta < nu
    c = 0;   % do NOT transmit
  else
    c = 1;   % transmit
  end

end