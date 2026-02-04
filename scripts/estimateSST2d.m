%% Simulations of the 2D Single-Target Tracking Plant

clear;
close all;
rng(42);

steps = 100;
x0 = [1 0]';

plant = singleTargetTrackingPlant(Ts);

% Initialize vectors for storing states and outputs (include t=0).
X = zeros(plant.n, steps + 1);
Y = zeros(plant.p, steps + 1);

X(:,1) = x0;
Y(:,1) = plant.outputEq(x0);
for t = 2:steps + 1
  [Y(:,t), X(:,t)] = plant.update(X(:, t-1));
end

t = (0:steps) * Ts;
figure
plot(t, X');
title("Simulated State of the Actual Model")
xlabel('Time (s)');

figure
subplot(2, 1, 1)
plot(t, Y(1:2,:)');
title("Simulated Output of nodes 1 and 2")
subplot(2, 1, 2)
plot(t, Y(3:4,:)');
title("Simulated Output of nodes 3 and 4")
xlabel('Time (s)');
