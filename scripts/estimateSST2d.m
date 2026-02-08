%% Simulations of the 2D Single-Target Tracking Plant

clear;
close all;
rng(42);

%% Parameters

T = 1000;  % Number of Simulation Steps
Ts = 0.1;  % Sampling Period
outputNoiseStd = 10;

x0 = [14 14 1800 2000]';  % v_x, v_y, p_x, p_y

nodeCount = 100;
sensorCount = 20;
maxLength = 5000;

%% Network Definition

[netGraph] = createSpatialNetwork(nodeCount, sensorCount, maxLength);

%% Model Simulation

plant = SingleTarget2dModel(Ts, sensorCount, outputNoiseStd, T);
plant.simulate(x0);

%% Estimators

% TODO: Review initialization
x0_hat = x0;
P0 = diag([1e2 1e2 1e6 1e6]);

ckf = CKF(plant, Ts, T);
ckf.run(x0_hat, P0);

%% Plotting

plotNetwork(netGraph, maxLength);

plant.plotTrajectory();
plant.plotOutputs();

ckf.plotTrajectory();
ckf.plotRmse();
