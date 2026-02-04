%% Description
% This script simulates the development of a system where n_nodes sensors
% monitor a system a communicate information to eachother based on a
% deterministic criteria. Each of the sensors get one row in the C matrix as
% their output.
% The Kalman filter is implemented in the information form to ease
% transition to a robust Kalman filter (minimax filter)
% The plant is a standard state-space model.
% All nodes are assumed to be connected.

clear;
close all;
rng(42);

%% Initialization

% Number of timesteps in simulation
simlength = 100;

% Number of nodes and states
n_nodes = 3;

n_states = 2;

% Initialize vectors for storing states, outputs etc.
Y = cell(simlength, 1); % True output from each node
X = cell(simlength, 1); % True states

% Estimated state as estimated by each node
X_hat = cell(simlength, n_nodes);

% Actual Plant State Space Model
A = [0 1; -0.8 -1];
B = [0.1; 0.2];
C = [2 0; 1 0; 0 1];
D = [0.5; 0.1; 0.1];
plant = StateSpaceModel(A, B, C, D);

x0 = [5; 0]; % Plant Initial Conditions

% Kalman Filter parameters
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

% For starters assuming all nodes are connected
pi = ones(n_nodes) / (n_nodes);

% The uncertainty factor, which limits the influence of error in stored estimates,
% when no transmission is made
delta = 0.05;

% State estimation error tolerance for sending information
epsilon = 0.05;
% Transmitting variable (1 if transmitting 0 if not)
c = zeros(n_nodes, 1);

%% Running the algorithm

% Iterate for all time steps
for t = 1:simlength
  % disp('new iteration')
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
    [q_bar, Omega_bar] = idkfPredict(q_fused, Omega_fused, A_hat, Q);

    % Correction step using information form Kalman Filter
    [q{t, i, i}, Omega{t, i, i}] = idkfCorrect(q_bar, Omega_bar, C_hat(i, :), R(i), Y{t}(i));

    % Return from information variables to states
    [X_hat{t, i}, x_hat_bar] = idkfGetState(Omega{t, i, i}, q{t, i, i}, Omega_bar, q_bar);

    % Evaluate whether to send data to other sensors for next iteration:
    if t == 1
      c(i) = 1;
    else
      c(i) = checkDeterministicFusionConditions(X_hat{t, i}, x_hat_bar, Omega{t, i, i}, epsilon);
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
