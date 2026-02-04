%% Check if a node should transmit its data to the other nodes
function c = checkStochasticFusionConditions(X, x_bar, Z)
  % Calculate the error in state estimates
  error = X - x_bar;

  % Implement eq. 10 in stochastic-event-triggered, but with states in
  % stead of outputs
  nu = exp(-1/2 * error' * Z * error);

  % Generate random number
  zeta = rand;

  % If error is larger than tolerance, then transmit
  if zeta < nu
    c = 0;
  else
    c = 1;
  end
end
