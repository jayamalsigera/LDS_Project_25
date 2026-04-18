
%% Simulations of the 2D Single-Target Tracking Plant

clear;
close all;
clc;
rng(42);


%% Parameters

% T = 2500; % Number of Simulation Steps
T = 1000;
% T = 100;
% T = 10;
% T = 2;

Ts = 0.1; % Sampling Period
outputNoiseStd = 10;

x0 = [25 25 1000 1000]'; % Initial true state: [vx vy px py]
% nodeCount = 100;% Total nodes in the network
% sensorCount = 20; % Number of sensor nodes
nodeCount = 20;
sensorCount = 4;

maxLength = 5000; % Spatial extent for random network generation

consensusSteps = 3; % Total Iterations for DSEACP consensus (L)
% consensusSteps = 50;

% Event-trigger parameters for DKF/RDKF/SRDKF
dkfAlpha = 10;
% dkfAlpha = 1;
dkfBeta = 0.2;
dkfDelta = 0.5;

% KL-divergence tolerance for robust filters (b)
klTolerance = 0.05;

% Stochastic Trigger error norm weights Matrix (Z)
errorNormWeightsClosed = 300 * eye(length(x0));
errorNormWeightsOpen = 300 * eye(2);
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
x0_hat = x0;  %Initial estimate
 P0 = 1e3 * eye(size(x0, 1)); % No Prior
%P0 = 1e-6 * eye(size(x0, 1));  % "Perfect Knowledge"

% Initiate all filters
ckf = CKF(plant, Ts, T);
dseacp = DSEACP(plant, Ts, T, netGraph, consensusSteps);
dkf = DKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta);
rdkf = RDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, klTolerance);
srdkfClosed = SRDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                    klTolerance, errorNormWeightsClosed, 'closed');

srdkfOpen = SRDKF(plant, Ts, T, netGraph, dkfAlpha, dkfBeta, dkfDelta, ...
                  klTolerance, errorNormWeightsOpen, 'open');
%% Monte Carlo simulations

% totalRuns = 200;
% totalRuns = 100;
totalRuns = 10;
 %totalRuns = 2;

% Preallocate RMSE and Transmission rate logs
ckfRmse = zeros(totalRuns, T + 1);
dseacpRmse = zeros(totalRuns, T + 1);
dkfRmse = zeros(totalRuns, T + 1);
rdkfRmse = zeros(totalRuns, T + 1);
srdkfClosedRmse = zeros(totalRuns, T + 1);
srdkfOpenRmse   = zeros(totalRuns, T + 1);


dkfTxRate = zeros(totalRuns, T + 1);
rdkfTxRate = zeros(totalRuns, T + 1);
srdkfClosedTxRate = zeros(totalRuns, T + 1);
srdkfOpenTxRate   = zeros(totalRuns, T + 1);

h = waitbar(0, 'Running simulations');
tic
for run = 1:totalRuns
  % Simulate one trajectory of the plant
  mdlSample = plant.simulate(x0);

  % Run each estimator on the same data and log RMSE + transmission rate
  ckfSample = ckf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
  ckfRmse(run, :) = ckfSample.RMSE;

  dseacpSample = dseacp.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
  dseacpRmse(run, :) = dseacpSample.RMSE;

  dkfSample = dkf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
  dkfRmse(run, :) = dkfSample.RMSE;
  dkfTxRate(run, :) = dkfSample.txRate;

  rdkfSample = rdkf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
  rdkfRmse(run, :) = rdkfSample.RMSE;
  rdkfTxRate(run, :) = rdkfSample.txRate;

  srdkfClosedSample = srdkfClosed.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
  srdkfClosedRmse(run, :) = srdkfClosedSample.RMSE;
  srdkfClosedTxRate(run, :) = srdkfClosedSample.txRate;

  srdkfOpenSample = srdkfOpen.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
  srdkfOpenRmse(run, :) = srdkfOpenSample.RMSE;
  srdkfOpenTxRate(run, :) = srdkfOpenSample.txRate;

  waitbar(run / totalRuns, h, sprintf('Run %d/%d', run, totalRuns));
end
fprintf('Elapsed: %.3f s\n', toc);
close(h)

%% Plotting

if true
  % if false
  plotNetwork(netGraph, maxLength); % Visualize node layout

  mdlSample.plotTrajectory(); % True trajectory
  mdlSample.plotOutputs(); % Sensor outputs

  % Estimated trajectories
  ckfSample.plotTrajectory(mdlSample.X);
  dseacpSample.plotTrajectory(mdlSample.X);
  dkfSample.plotTrajectory(mdlSample.X);
  rdkfSample.plotTrajectory(mdlSample.X);
  srdkfClosedSample.plotTrajectory(mdlSample.X);
  srdkfOpenSample.plotTrajectory(mdlSample.X);

  % RMSE comparison
  figure
  t = (0:T) * Ts;
  semilogy(t, mean(ckfRmse, 1), 'DisplayName', 'CKF');
  hold on;
  semilogy(t, mean(dseacpRmse, 1), 'DisplayName', 'DSEA-CP (L=3)');
  semilogy(t, mean(dkfRmse, 1), 'DisplayName', 'DKF');
  semilogy(t, mean(rdkfRmse, 1), 'DisplayName', 'RDKF');
  semilogy(t, mean(srdkfClosedRmse, 1), 'DisplayName', 'SRDKF-Closed');
  semilogy(t, mean(srdkfOpenRmse, 1), 'DisplayName', 'SRDKF-Open');
  hold off;
  title("RMSE vs Time");
  xlabel('Time (s)');
  ylabel('RMSE');
  legend();
  grid();


  % Transmission rate comparison
  figure
  t = (0:T) * Ts;
  plot(t, mean(dkfTxRate, 1), 'DisplayName', 'DKF');
  hold on
  plot(t, mean(rdkfTxRate, 1), 'DisplayName', 'RDKF');
  plot(t, mean(srdkfClosedTxRate, 1), 'DisplayName', 'SRDKF-Closed');
  plot(t, mean(srdkfOpenTxRate, 1), 'DisplayName', 'SRDKF-Open');
  hold off
  title("TX Rate vs Time");
  xlabel('Time (s)');
  ylabel('TX Rate');
  legend();
  grid();
end
