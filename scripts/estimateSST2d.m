%% Simulations of the 2D Single-Target Tracking Plant

clear;
close all;
rng(42);

%% Parameters

% T = 1000; % Number of Simulation Steps
T = 200; % Number of Simulation Steps
% T = 10; % Number of Simulation Steps
Ts = 0.1; % Sampling Period
outputNoiseStd = 10;

x0 = [50 50 1800 2000]'; % v_x, v_y, p_x, p_y

% nodeCount = 100;
% sensorCount = 20;
nodeCount = 20;
sensorCount = 4;

maxLength = 5000;

consensusSteps = 3;
% consensusSteps = 50;

dkfAlpha = 10;
dkfBeta = 0.2;
dkfDelta = 0.5;

%% Network Definition

netGraph = createSpatialNetwork(nodeCount, sensorCount, maxLength);

% TODO: Connectivity check
% TODO: A doubly stochastic check

%% Model Simulation

plant = SingleTarget2dModel(Ts, sensorCount, outputNoiseStd, T);

% TODO: Stabilizability check
% TODO: Detectability (or Observability?) check

%% Estimators

% TODO: Review initialization
x0_hat = x0;
% P0 = 1e12 * eye(size(x0, 1));  % No Prior
P0 = 1e-12 * eye(size(x0, 1));  % "Perfect Knowledge"

ckf = CKF(plant, Ts, T);
dseacp = DSEACP(plant, Ts, T, netGraph, consensusSteps);
dkf = DKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta);

%% Monte Carlo simulations

% totalRuns = 200;
totalRuns = 100;
% totalRuns = 10;
% totalRuns = 1;

ckfRmse = zeros(totalRuns, T + 1);
dseacpRmse = zeros(totalRuns, T + 1);
dkfRmse = zeros(totalRuns, T + 1);

h = waitbar(0, 'Running simulations');
for run = 1:totalRuns
  mdlSample = plant.simulate(x0);

  ckfSample = ckf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
  ckfRmse(run, :) = ckfSample.RMSE;

  % dseacpSample = dseacp.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
  % dseacpRmse(run, :) = dseacpSample.RMSE;

  dkfSample = dkf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
  dkfRmse(run, :) = dkfSample.RMSE;

  waitbar(run / totalRuns, h, sprintf('Run %d/%d', run, totalRuns));
end
close(h)

%% Plotting

if true
% if false
  plotNetwork(netGraph, maxLength);

  mdlSample.plotTrajectory();
  mdlSample.plotOutputs();

  ckfSample.plotTrajectory(mdlSample.X);
  dseacpSample.plotTrajectory(mdlSample.X);
  dkfSample.plotTrajectory(mdlSample.X);

  figure
  t = (0:T) * Ts;
  semilogy(t, mean(ckfRmse, 1), 'DisplayName', 'CKF');
  hold on;
  % semilogy(t, mean(dseacpRmse, 1), 'DisplayName', 'DSEA-CP (L=3)');
  semilogy(t, mean(dkfRmse, 1), 'DisplayName', 'DKF');
  hold off;
  title("RMSE vs Time");
  xlabel('Time (s)');
  ylabel('RMSE');
  legend();
  grid();
end
