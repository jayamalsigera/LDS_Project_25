%% Simulations of the 2D Single-Target Tracking Plant

clear;
close all;
rng(42);

%% Parameters

T = 1000;  % Number of Simulation Steps
Ts = 0.1;  % Sampling Period

x0 = [14 14 1800 2000]';  % v_x, v_y, p_x, p_y

nodeCount = 100;
sensorCount = 20;
maxLength = 5000;

%% Network Definition

[netGraph] = createSpatialNetwork(nodeCount, sensorCount, maxLength);

%% Model Simulation

plant = singleTarget2dModel(Ts, sensorCount, T);
plant.simulate(x0);

%% Estimators

% TODO: Review initialization
x0_hat = x0;
P0 = diag([1e2 1e2 1e6 1e6]);

ckf = CKF(plant, T);
ckf.run(x0_hat, P0);

%% Plotting

plotNetwork(netGraph, maxLength);

figure
plot(plant.X(3, :), plant.X(4, :));
title("Simulated Trajectory (State)")
xlabel('$p_x$', 'Interpreter', 'latex');
ylabel('$p_y$', 'Interpreter', 'latex');
grid()
xlim([0, maxLength]);
ylim([0, maxLength]);

figure
plot(ckf.x_hat(3, :), ckf.x_hat(4, :));
title("CKF Estimated Trajectory")
xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
grid()
xlim([0, maxLength]);
ylim([0, maxLength]);

figure
t = (0:T) * Ts;
subplot(2, 1, 1)
plot(t, plant.Y(1:2,:)');
title("Simulated Output of Sensor 1")
subplot(2, 1, 2)
plot(t, plant.Y(39:40,:)');
title("Simulated Output of Sensor 20")
xlabel('Time (s)');
