
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
disp("Creating Network")
netGraph = createSpatialNetwork(nodeCount, sensorCount, maxLength);

assertConnected(netGraph);
% TODO: A doubly stochastic check

%% Model Simulation
disp("Simulating target dynamics")
plant = SingleTarget2dModel(Ts, sensorCount, outputNoiseStd, T);

% TODO: Stabilizability check
% TODO: Detectability (or Observability?) check

%% Estimators
disp("Initializing estimators")

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
disp("Running Monte Carlo simulations")

% totalRuns = 200;
% totalRuns = 100;
totalRuns = 10;
% totalRuns = 2;

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

tic

% Showcase run (serial) — keeps sample objects available for trajectory plots
mdlSample = plant.simulate(x0);
ckfSample = ckf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
dseacpSample = dseacp.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
dkfSample = dkf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
rdkfSample = rdkf.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
srdkfClosedSample = srdkfClosed.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);
srdkfOpenSample = srdkfOpen.estimate(x0_hat, P0, mdlSample.X, mdlSample.Y);

ckfRmse(1, :) = ckfSample.RMSE;
dseacpRmse(1, :) = dseacpSample.RMSE;
dkfRmse(1, :) = dkfSample.RMSE;
dkfTxRate(1, :) = dkfSample.txRate;
rdkfRmse(1, :) = rdkfSample.RMSE;
rdkfTxRate(1, :) = rdkfSample.txRate;
srdkfClosedRmse(1, :) = srdkfClosedSample.RMSE;
srdkfClosedTxRate(1, :) = srdkfClosedSample.txRate;
srdkfOpenRmse(1, :) = srdkfOpenSample.RMSE;
srdkfOpenTxRate(1, :) = srdkfOpenSample.txRate;

q = parforWaitbar(totalRuns - 1, 'Running simulations');
parfor run = 2:totalRuns
  % Simulate one trajectory of the plant
  s = plant.simulate(x0);

  % Run each estimator on the same data and log RMSE + transmission rate
  ckfRmse(run, :) = ckf.estimate(x0_hat, P0, s.X, s.Y).RMSE;
  dseacpRmse(run, :) = dseacp.estimate(x0_hat, P0, s.X, s.Y).RMSE;

  dkfRun = dkf.estimate(x0_hat, P0, s.X, s.Y);
  dkfRmse(run, :) = dkfRun.RMSE;
  dkfTxRate(run, :) = dkfRun.txRate;

  rdkfRun = rdkf.estimate(x0_hat, P0, s.X, s.Y);
  rdkfRmse(run, :) = rdkfRun.RMSE;
  rdkfTxRate(run, :) = rdkfRun.txRate;

  srdkfClosedRun = srdkfClosed.estimate(x0_hat, P0, s.X, s.Y);
  srdkfClosedRmse(run, :) = srdkfClosedRun.RMSE;
  srdkfClosedTxRate(run, :) = srdkfClosedRun.txRate;

  srdkfOpenRun = srdkfOpen.estimate(x0_hat, P0, s.X, s.Y);
  srdkfOpenRmse(run, :) = srdkfOpenRun.RMSE;
  srdkfOpenTxRate(run, :) = srdkfOpenRun.txRate;

  send(q, run);
end
fprintf('Elapsed: %.3f s\n', toc);

%% Plotting

plottingEnabled = true;
if plottingEnabled
  disp("Plotting results.")
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

  % Consistent color per estimator across all plots
  colors = struct( ...
    'CKF',          [0.00 0.45 0.74], ...
    'DSEACP',       [0.85 0.33 0.10], ...
    'DKF',          [0.93 0.69 0.13], ...
    'RDKF',         [0.49 0.18 0.56], ...
    'SRDKFClosed',  [0.47 0.67 0.19], ...
    'SRDKFOpen',    [0.30 0.75 0.93]);

  % RMSE comparison
  figure
  t = (0:T) * Ts;
  semilogy(t, mean(ckfRmse, 1), 'Color', colors.CKF, 'DisplayName', 'CKF');
  hold on;
  semilogy(t, mean(dseacpRmse, 1), 'Color', colors.DSEACP, 'DisplayName', 'DSEA-CP (L=3)');
  semilogy(t, mean(dkfRmse, 1), 'Color', colors.DKF, 'DisplayName', 'DKF');
  semilogy(t, mean(rdkfRmse, 1), 'Color', colors.RDKF, 'DisplayName', 'RDKF');
  semilogy(t, mean(srdkfClosedRmse, 1), 'Color', colors.SRDKFClosed, 'DisplayName', 'SRDKF-Closed');
  semilogy(t, mean(srdkfOpenRmse, 1), 'Color', colors.SRDKFOpen, 'DisplayName', 'SRDKF-Open');
  hold off;
  title("RMSE vs Time");
  xlabel('Time (s)');
  ylabel('RMSE');
  legend();
  grid();


  % Transmission rate comparison
  figure
  t = (0:T) * Ts;
  plot(t, mean(dkfTxRate, 1), 'Color', colors.DKF, 'DisplayName', 'DKF');
  hold on
  plot(t, mean(rdkfTxRate, 1), 'Color', colors.RDKF, 'DisplayName', 'RDKF');
  plot(t, mean(srdkfClosedTxRate, 1), 'Color', colors.SRDKFClosed, 'DisplayName', 'SRDKF-Closed');
  plot(t, mean(srdkfOpenTxRate, 1), 'Color', colors.SRDKFOpen, 'DisplayName', 'SRDKF-Open');
  hold off
  title("TX Rate vs Time");
  xlabel('Time (s)');
  ylabel('TX Rate');
  legend();
  grid();
end
