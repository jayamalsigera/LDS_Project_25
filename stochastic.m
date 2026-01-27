%% Description
% This script simulates the development of a system where n_nodes sensors
% monitor a system a communicate information to each other based on a
% stochastic criteria.
%
% Each of the sensors get one row in the C matrix as their individual C
%
% The Kalman filter is implemented in the information form to ease
% transition to a robust Kalman filter (minimax filter)
%
% The plant is a standard state-space model.

clear;
close all;
rng(42);

%% Initialization

% Number of timesteps in simulation
simlength = 10;

% Number of nodes and states
n_nodes = 3;

n_states = 2;

% Initialize vectors for storing states, outputs etc.
Y = cell(simlength, 1); % True output from each node
X = cell(simlength, 1); % True states

% Estimated state as estimated by each node
X_hat = cell(simlength, n_nodes);

% TODO: What are these doing?
c_test = X_hat;
nu = c_test;
zeta = nu;
jimmy = nu;
c_tune = nu;

% Actual Plant State Space Model
A = [0 1; -0.8 -1];
B = [0.1; 0.2];
C = [2 0; 1 0; 0 1];
D = [0.5; 0.1; 0.1];
plant = StateSpaceModel(A, B, C, D);

x0 = [5; 0]; % Plant Initial Conditions

% Kalman Filter Parameters
A_hat = A;
C_hat = C;
Q = 0.1 * eye(n_states);
R = D .* D;
x0_hat = x0; % Filter Initial Conditions

% Variables for information form Kalman Filter
% The vector q interpreted as for everytime t the vector q as seen by each state as
% predicted by each state (q(t,2,1) is what node 1 thinks node 2's q is)
% The spot (t,i,i) contains the nodes knowledge of its own information
% vector and is always correct
q = cell(simlength, n_nodes, n_nodes);

% The same convention is used in Omega
Omega = cell(simlength, n_nodes, n_nodes);

for i = 1:n_nodes
  % Omega is initialized as the identity matrix
  Omega{1, i, i} = eye(n_states);
  q{1, i, i} = Omega{1, i, i} \ x0_hat;
end

for i = 1:n_nodes
  for j = 1:n_nodes
    if i ~= j
      Omega{1, i, j} = eye(n_states); % or some prior
      q{1, i, j} = zeros(n_states, 1); % or x0_hat, or q{1,j,j}
    end
  end
end

% Sending and receiving information

% The connections between different nodes are defined in the matrix below
% The diagonal must be 1's
connections = [1 1 0
               1 1 1
               0 1 1];

% Define pi, which contains the weights for fusing data
pi = zeros(size(connections));

% Each element in pi is given the weight: number of connections + 1
% Since a node is also counted as connected to itself this is the sum of
% elements in a row
for i = 1:length(connections(:, 1))
  pi(i, :) = connections(i, :) / sum(connections(i, :));
end

% The uncertainty factor, which limits the influence of error in stored estimates,
% when no transmission is made
delta = 0.05;

% The matrix used for calculting the exponential in stochastic
% event-triggered schedule
Z = eye(n_states);

% State estimation error tolerance for sending information
epsilon = 0.05;

% Transmitting variable (1 if transmitting 0 if not)
c = zeros(n_nodes, 1);

%% Running the Algorithm

% Iterate for all time steps
for t = 1:simlength
  % Update real state:
  if t == 1
    [Y{t}, X{t}] = plant.update(x0);
  else
    [Y{t}, X{t}] = plant.update(X{t - 1});
  end

  % Can definitely be done more tidy but this seems a nice way to start
  % where the time t is intuitive

  % For each node
  for i = 1:n_nodes
    % Start by fusing the information
    if t == 1
      % For first iteration take the initialized variables, which
      % contain an initial guess
      q_fused = q{t, i, i};
      Omega_fused = Omega{t, i, i};

    else
      % For all subsequent steps the fusing of data is done using the
      % fuseData Function
      [q_new, Omega_new, q_fused, Omega_fused] = fuseData(reshape(q(t - 1, :, :), n_nodes, n_nodes), reshape(Omega(t - 1, :, :), n_nodes, n_nodes), c_prev, pi, i, delta);
      q(t, :, i) = q_new.';
      Omega(t, :, i) = Omega_new.';
    end

    % Prediction step using information form Kalman Filter
    [q_bar, Omega_bar] = kfPredict(q_fused, Omega_fused, A_hat, Q);

    % Correction step using information form Kalman Filter
    [q{t, i, i}, Omega{t, i, i}] = kfCorrect(q_bar, Omega_bar, C_hat(i, :), R(i), Y{t}(i));

    % Return from information variables to states
    [X_hat{t, i}, x_hat_bar] = getStates(Omega{t, i, i}, q{t, i, i}, Omega_bar, q_bar);

    % Evaluate whether to send data to other sensors for next iteration:
    if t == 1
      c(i) = 1;
    else
      % c(i) =checkFusionCondition(X_hat{t,i},x_hat_bar,Omega{t,i,i},epsilon);
      c(i) = checkFusionConditionStochastic(X_hat{t, i}, x_hat_bar, eye(n_states));

    end

  end

  c_prev = c;
end

%% Plotting of the resulting state estimates

figure(1)
% Returning from cell array to actual array for plotting
X_plot = zeros(simlength, n_states);

for t = 1:simlength
  X_plot(t, :) = X{t}.';
end

X_hat_plot = zeros(simlength, n_states, n_nodes);

for t = 1:simlength

  for i = 1:n_nodes
    X_hat_plot(t, :, i) = X_hat{t, i}.';
  end

end

% Plotting of actual and estimated states for each sensor:
plot(X_plot(:, 1))
hold on
plot(X_hat_plot(:, 1, 1))
plot(X_hat_plot(:, 1, 2))
plot(X_hat_plot(:, 1, 3))
legend('True X', 'Node 1', 'Node 2', 'Node 3')
hold off
figure(2)
plot(X_plot(:, 2))
hold on
plot(X_hat_plot(:, 2, 1))
plot(X_hat_plot(:, 2, 2))
plot(X_hat_plot(:, 2, 3))
legend('True X', 'Node 1', 'Node 2', 'Node 3')
hold off

%% Functions

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


% Returning from information form to statespace-form state vectors
function [X, x_bar] = getStates(Omega, q, Omega_bar, q_bar)
  % Modification of first equation in III B in:
  % https://isas.iar.kit.edu/pdf/Fusion17_Pfaff-IDKF.pdf
  X = Omega \ q;
  x_bar = Omega_bar \ q_bar;
end

% Check if a node should transmit its data to the other nodes
function c = checkFusionCondition(X, x_bar, Omega, epsilon)
  error = X - x_bar;

  % If error is larger than tolerance, then transmit
  if error' * (Omega \ error) < epsilon
    c = 0;
  else
    c = 1;
  end

end

% Check if a node should transmit its data to the other nodes
function c = checkFusionConditionStochastic(X, x_bar, Z)
  % Calculate the erro in state estimates
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
