%% Check if a node should transmit its data to the other nodes
function c = checkDeterministicFusionConditions(X, x_bar, Omega, epsilon)
  error = X - x_bar;

  % If error is larger than tolerance, then transmit
  if error' * (Omega \ error) < epsilon
    c = 0;
  else
    c = 1;
  end

end
