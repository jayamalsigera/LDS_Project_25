%% Simulations of the 2D Single-Target Tracking Plant

clear;
close all;
rng(42);

%% Parameters

steps = 1000;
Ts = 0.1;

x0 = [14 14 1800 2000]';  % v_x, v_y, p_x, p_y

nodeCount = 100;
sensorCount = 20;
maxLength = 5000;

%% Plant Definition

[netGraph] = createSpatialNetwork(nodeCount, sensorCount, maxLength);
plant = singleTargetTracking2dPlant(Ts, sensorCount);

%% Simulation Loop

% Initialize vectors for storing states and outputs (include t=0).
X = zeros(plant.n, steps + 1);
Y = zeros(plant.p, steps + 1);

X(:,1) = x0;
Y(:,1) = plant.outputEq(x0);
for t = 2:steps + 1
  [Y(:,t), X(:,t)] = plant.update(X(:, t-1));
end

%% Plotting

plotNetwork(netGraph, maxLength);

figure
plot(X(3, :), X(4, :));
title("Simulated Trajectory (State)")
xlabel('$p_x$', 'Interpreter', 'latex');
ylabel('$p_y$', 'Interpreter', 'latex');
grid()
xlim([0, maxLength]);
ylim([0, maxLength]);

figure
t = (0:steps) * Ts;
subplot(2, 1, 1)
plot(t, Y(1:2,:)');
title("Simulated Output of Sensor 1")
subplot(2, 1, 2)
plot(t, Y(39:40,:)');
title("Simulated Output of Sensor 20")
xlabel('Time (s)');
