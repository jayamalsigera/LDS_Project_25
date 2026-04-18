%% Simulations of the 2D Single-Target Tracking Plant under the Least-Favorable Model
%
% Mirrors estimateSST2d.m, except the per-run trajectories (X, Y) are drawn from
% the least-favorable model (Levy & Nikoukhah 2013, Section V) at tolerance b,
% rather than from the nominal plant. This is the setup used by Ghion & Zorzi
% (2023) in their Monte Carlo study.

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

turnRate = 0; % Straight Line
% turnRate = 0.005; % Coordinated-turn rate [rad/s]; set to 0 for straight line

x0 = [25 25 1000 1000]'; % Initial true state: [vx vy px py]
% nodeCount = 100;% Total nodes in the network
% sensorCount = 20; % Number of sensor nodes
nodeCount = 20;
sensorCount = 4;

maxLength = 5000; % Spatial extent for random network generation

consensusSteps = 3; % Total Iterations for DSEACP consensus (L)
% consensusSteps = 50;

% Event-trigger parameters for DKF/RDKF/SRDKF
% dkfAlpha = 1;
dkfAlpha = 10;
% dkfBeta = 0.2;
dkfBeta = 1;
% dkfDelta = 0.1;
dkfDelta = 1;

% KL-divergence tolerance for robust filters (b). The LFM uses the same b.
klTolerance = 0.1;

% Stochastic Trigger error norm weights Matrix (Z)
errorNormWeightsClosed = 300 * eye(length(x0));
errorNormWeightsOpen = 300 * eye(2);
%% Network Definition
disp("Creating Network")
netGraph = createSpatialNetwork(nodeCount, sensorCount, maxLength);

assertConnected(netGraph);

%% Model Simulation
disp("Simulating target dynamics")
plant = SingleTarget2dModel(Ts, sensorCount, outputNoiseStd, T, turnRate);

assertStabilizable(plant.A, plant.B);
assertDetectable(plant.A, plant.C);

%% Estimators
disp("Initializing estimators")

% TODO: Review initialization
x0_hat = x0;  %Initial estimate
 P0 = 1e3 * eye(size(x0, 1)); % No Prior
%P0 = 1e-6 * eye(size(x0, 1));  % "Perfect Knowledge"

% Least-favorable data generator. Forward/backward sweeps (G_t, H_t, L_t) are
% computed once here; only the per-run randn draws happen inside the parfor.
disp("Precomputing least-favorable model")
lfm = LeastFavorableModel(plant, P0, klTolerance, T);

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
disp("Running Monte Carlo simulations (LFM data)")

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

% Showcase run (serial) — keeps sample objects available for trajectory plots.
% Also draw a nominal-model sample for side-by-side visual comparison.
mdlSample = lfm.simulate(x0);
nomSample = plant.simulate(x0);
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
  % Draw one trajectory from the least-favorable model
  s = lfm.simulate(x0);

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

  % True trajectory (LFM sample) overlaid with a nominal-model sample, so
  % the perturbation introduced by the least-favorable model is visible.
  figure
  plot(nomSample.X(3, :), nomSample.X(4, :), 'DisplayName', 'Nominal');
  hold on
  plot(mdlSample.X(3, :), mdlSample.X(4, :), 'DisplayName', 'LFM');
  hold off
  title("Simulated Trajectory")
  xlabel('$p_x$', 'Interpreter', 'latex');
  ylabel('$p_y$', 'Interpreter', 'latex');
  legend();
  grid()

  % Estimated trajectories (estimator mean vs LFM truth; nominal overlaid
  % as a visual reference for how far the worst-case model drifts).
  plotEstimatedVsTruth = @(est, label) plotWithNominal( ...
    est, label, mdlSample.X, nomSample.X);
  plotEstimatedVsTruth(ckfSample,         'CKF');
  plotEstimatedVsTruth(dseacpSample,      'DSEA-CP');
  plotEstimatedVsTruth(dkfSample,         'DKF');
  plotEstimatedVsTruth(rdkfSample,        'RDKF');
  plotEstimatedVsTruth(srdkfClosedSample, 'SRDKF-Closed');
  plotEstimatedVsTruth(srdkfOpenSample,   'SRDKF-Open');

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
  title("RMSE vs Time (LFM data)");
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
  title("TX Rate vs Time (LFM data)");
  xlabel('Time (s)');
  ylabel('TX Rate');
  legend();
  grid();
end

function plotWithNominal(estSample, label, Xlfm, Xnom)
  if isprop(estSample, 'X_hat')
    % Distributed estimators: n x N x (T+1). Collapse node dim.
    meanX_hat = squeeze(mean(estSample.X_hat, 2));
  else
    % Centralized estimators (e.g., CKF): n x (T+1), no node dim.
    meanX_hat = estSample.x_hat;
  end
  figure
  plot(meanX_hat(3, :), meanX_hat(4, :), 'DisplayName', label);
  hold on
  plot(Xlfm(3, :), Xlfm(4, :), 'DisplayName', 'LFM truth');
  plot(Xnom(3, :), Xnom(4, :), '--', 'DisplayName', 'Nominal');
  hold off
  title(sprintf("%s Estimated Trajectory", label));
  xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
  ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
  legend();
  grid();
end
