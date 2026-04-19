function plotSST2dRun(path)
%PLOTSST2DRUN  Re-create figures for a saved estimateSST2d[Lfm] run.
%
% Works for both the nominal-plant script (estimateSST2d) and the
% least-favorable model script (estimateSST2dLfm); the LFM branch is
% taken when runData.samples has a nomSample field.

  runData = loadRun(path);

  params   = runData.params;
  netGraph = runData.netGraph;
  results  = runData.results;
  samples  = runData.samples;

  T  = params.T;
  Ts = params.Ts;
  isLfm = isfield(samples, 'nomSample');

  % Network layout
  plotNetwork(netGraph, params.maxLength);

  % Truth trajectory
  if isLfm
    figure
    plot(samples.nomSample.X(3, :), samples.nomSample.X(4, :), ...
         'DisplayName', 'Nominal');
    hold on
    plot(samples.mdlSample.X(3, :), samples.mdlSample.X(4, :), ...
         'DisplayName', 'LFM');
    hold off
    title("Simulated Trajectory")
    xlabel('$p_x$', 'Interpreter', 'latex');
    ylabel('$p_y$', 'Interpreter', 'latex');
    legend(); grid()
  else
    samples.mdlSample.plotTrajectory();
    samples.mdlSample.plotOutputs();
  end

  % Estimated trajectories per filter
  estLabels = {'CKF',        'DSEA-CP',     'DKF',        'RDKF', ...
               'SRDKF-Closed',               'SRDKF-Open'};
  estFields = {'ckfSample',  'dseacpSample','dkfSample',  'rdkfSample', ...
               'srdkfClosedSample',          'srdkfOpenSample'};
  for k = 1:numel(estFields)
    if ~isfield(samples, estFields{k}); continue; end
    est = samples.(estFields{k});
    if isLfm
      plotEstVsTruthWithNominal(est, estLabels{k}, ...
                                samples.mdlSample.X, samples.nomSample.X);
    else
      est.plotTrajectory(samples.mdlSample.X);
    end
  end

  % Consistent color per estimator across all plots
  colors = struct( ...
    'CKF',          [0.00 0.45 0.74], ...
    'DSEACP',       [0.85 0.33 0.10], ...
    'DKF',          [0.93 0.69 0.13], ...
    'RDKF',         [0.49 0.18 0.56], ...
    'SRDKFClosed',  [0.47 0.67 0.19], ...
    'SRDKFOpen',    [0.30 0.75 0.93]);

  t = (0:T) * Ts;

  % RMSE comparison
  figure
  semilogy(t, mean(results.ckfRmse, 1),         'Color', colors.CKF,         'DisplayName', 'CKF');
  hold on;
  semilogy(t, mean(results.dseacpRmse, 1),      'Color', colors.DSEACP,      'DisplayName', 'DSEA-CP (L=3)');
  semilogy(t, mean(results.dkfRmse, 1),         'Color', colors.DKF,         'DisplayName', 'DKF');
  semilogy(t, mean(results.rdkfRmse, 1),        'Color', colors.RDKF,        'DisplayName', 'RDKF');
  semilogy(t, mean(results.srdkfClosedRmse, 1), 'Color', colors.SRDKFClosed, 'DisplayName', 'SRDKF-Closed');
  semilogy(t, mean(results.srdkfOpenRmse, 1),   'Color', colors.SRDKFOpen,   'DisplayName', 'SRDKF-Open');
  hold off;
  if isLfm
    title("RMSE vs Time (LFM data)");
  else
    title("RMSE vs Time");
  end
  xlabel('Time (s)'); ylabel('RMSE');
  legend(); grid();

  % Transmission rate comparison
  figure
  plot(t, mean(results.dkfTxRate, 1),         'Color', colors.DKF,         'DisplayName', 'DKF');
  hold on
  plot(t, mean(results.rdkfTxRate, 1),        'Color', colors.RDKF,        'DisplayName', 'RDKF');
  plot(t, mean(results.srdkfClosedTxRate, 1), 'Color', colors.SRDKFClosed, 'DisplayName', 'SRDKF-Closed');
  plot(t, mean(results.srdkfOpenTxRate, 1),   'Color', colors.SRDKFOpen,   'DisplayName', 'SRDKF-Open');
  hold off
  if isLfm
    title("TX Rate vs Time (LFM data)");
  else
    title("TX Rate vs Time");
  end
  xlabel('Time (s)'); ylabel('TX Rate');
  legend(); grid();
end

function plotEstVsTruthWithNominal(estSample, label, Xlfm, Xnom)
  if isprop(estSample, 'X_hat')
    meanX_hat = squeeze(mean(estSample.X_hat, 2));
  else
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
  legend(); grid();
end
