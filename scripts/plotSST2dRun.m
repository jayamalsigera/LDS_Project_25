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
  estLabels = {'CKF', 'CRKF', 'DSEA-CP', 'DKF', 'RDKF', ...
               'SDKF-Closed', 'SDKF-Open', 'SRDKF-Closed', 'SRDKF-Open'};
  estFields = {'ckfSample', 'crkfSample', 'dseacpSample', 'dkfSample', 'rdkfSample', ...
               'sdkfClosedSample', 'sdkfOpenSample', 'srdkfClosedSample', 'srdkfOpenSample'};
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

  % Colors grouped by estimator family; line styles separate variants within a family.
  %   Centralized  → blue family    (CKF solid, CRKF dashed)
  %   Consensus    → red/orange     (DSEA-CP solid)
  %   DKF family   → amber          (DKF solid, SDKF-Closed dashed, SDKF-Open dotted)
  %   RDKF family  → purple/magenta (RDKF solid, SRDKF-Closed dashed, SRDKF-Open dotted)
  c_ckf          = [0.00 0.45 0.74];   % MATLAB default blue
  c_crkf         = [0.30 0.65 0.90];   % lighter blue  (distinct from CKF)
  c_dseacp       = [0.85 0.33 0.10];   % MATLAB default orange-red
  c_dkf          = [0.93 0.69 0.13];   % MATLAB default amber
  c_rdkf         = [0.49 0.18 0.56];   % MATLAB default purple
  c_sdkf         = [0.10 0.75 0.65];   % teal  (SDKF family)
  c_srdkf        = [0.47 0.67 0.19];   % green (SRDKF family)

  t = (0:T) * Ts;

  lw = 1.3;   % uniform line width

  % RMSE comparison
  figure
  hold on
  semilogy(t, mean(results.ckfRmse,    1), '-',  'Color', c_ckf,    'LineWidth', lw, 'DisplayName', 'CKF');
  if isfield(results, 'crkfRmse')
    semilogy(t, mean(results.crkfRmse, 1), '--', 'Color', c_crkf,   'LineWidth', lw, 'DisplayName', 'CRKF');
  end
  semilogy(t, mean(results.dseacpRmse, 1), '-',  'Color', c_dseacp, 'LineWidth', lw, 'DisplayName', 'DSEA-CP');
  semilogy(t, mean(results.dkfRmse,    1), '-',  'Color', c_dkf,    'LineWidth', lw, 'DisplayName', 'DKF');
  semilogy(t, mean(results.rdkfRmse,   1), '-',  'Color', c_rdkf,   'LineWidth', lw, 'DisplayName', 'RDKF');
  if isfield(results, 'sdkfClosedRmse')
    semilogy(t, mean(results.sdkfClosedRmse, 1), '--', 'Color', c_sdkf,  'LineWidth', lw, 'DisplayName', 'SDKF-Closed');
    semilogy(t, mean(results.sdkfOpenRmse,   1), ':',  'Color', c_sdkf,  'LineWidth', lw+0.4, 'DisplayName', 'SDKF-Open');
  end
  semilogy(t, mean(results.srdkfClosedRmse, 1), '--', 'Color', c_srdkf, 'LineWidth', lw, 'DisplayName', 'SRDKF-Closed');
  semilogy(t, mean(results.srdkfOpenRmse,   1), ':',  'Color', c_srdkf, 'LineWidth', lw+0.4, 'DisplayName', 'SRDKF-Open');
  hold off
  title(sprintf('RMSE vs Time%s', optstr(isLfm, ' (LFM data)', '')));
  xlabel('Time (s)'); ylabel('RMSE');
  legend('Location', 'northeast'); grid();

  % Transmission rate comparison
  figure
  hold on
  plot(t, mean(results.dkfTxRate,  1), '-',  'Color', c_dkf,   'LineWidth', lw, 'DisplayName', 'DKF');
  plot(t, mean(results.rdkfTxRate, 1), '-',  'Color', c_rdkf,  'LineWidth', lw, 'DisplayName', 'RDKF');
  if isfield(results, 'sdkfClosedTxRate')
    plot(t, mean(results.sdkfClosedTxRate, 1), '--', 'Color', c_sdkf,  'LineWidth', lw, 'DisplayName', 'SDKF-Closed');
    plot(t, mean(results.sdkfOpenTxRate,   1), ':',  'Color', c_sdkf,  'LineWidth', lw+0.4, 'DisplayName', 'SDKF-Open');
  end
  plot(t, mean(results.srdkfClosedTxRate, 1), '--', 'Color', c_srdkf, 'LineWidth', lw, 'DisplayName', 'SRDKF-Closed');
  plot(t, mean(results.srdkfOpenTxRate,   1), ':',  'Color', c_srdkf, 'LineWidth', lw+0.4, 'DisplayName', 'SRDKF-Open');
  hold off
  title(sprintf('TX Rate vs Time%s', optstr(isLfm, ' (LFM data)', '')));
  xlabel('Time (s)'); ylabel('TX Rate');
  legend('Location', 'northeast'); grid();
end

function s = optstr(cond, a, b)
  if cond; s = a; else; s = b; end
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
