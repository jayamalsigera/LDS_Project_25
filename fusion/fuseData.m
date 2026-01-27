% Function to fuse the data between nodes
function [q_new, Omega_new, q_fused, Omega_fused] = fuseData(q, Omega, c, pi, i, delta)
  % We wish to create the summation eq. 10-11 from robust-kalman-eventtriggered

  % First data is copied from 1 object to the other, as the new objects
  % contain the same data as the old ones, unless updated
  Omega_new = Omega(:, i);
  q_new = q(:, i);

  % The lengths are redefined so as to not have to pass them
  n_nodes = length(c);
  n_states = length(q{i, i});
  % The sizes of the variables used to do the summation are defined
  sum_q = zeros(n_states, 1);
  sum_Omega = zeros(n_states);
  % Loop over all nodes except i itself
  for j = 1:n_nodes
    % skip the case where j=i
    if j == i
      continue
    end

    if c(j) == 1
      % If c = 1 take the new input
      q_new(j) = q(j, j);
      Omega_new(j) = Omega(j, j);
      % Here the new measurements are added to the summation object only gained by pi
      sum_q = sum_q + pi(i, j) * q{j, j};
      sum_Omega = sum_Omega + pi(i, j) * Omega{j, j};
    else
      % If c = 0 no transmission is received and stored estimate
      % is added scaled by pi and (1+delta)^-1 in stead of fresh
      % data
      sum_q = sum_q + 1 / (1 + delta) * pi(i, j) * q{j, i};
      sum_Omega = sum_Omega + 1 / (1 + delta) * pi(i, j) * Omega{j, i};

    end

  end

  % Information fusion
  q_fused = pi(i, i) * q{i, i} + sum_q;
  Omega_fused = pi(i, i) * Omega{i, i} + sum_Omega;
end
