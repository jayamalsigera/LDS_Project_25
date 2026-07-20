function plotSST2dTrajectories(path)
%PLOTSST2DTRAJECTORIES  Plot truth and estimated trajectories for a saved run.
%
% Works for both the nominal-plant script (estimateSST2d) and the
% least-favorable model script (estimateSST2dLfm); the LFM branch is
% taken when runData.samples has a nomSample field.

  runData = loadRun(path);

  samples = runData.samples;
  isLfm   = isfield(samples, 'nomSample');

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
               'SDKF-Closed', 'SDKF-Open', 'SRDKF-Closed', 'SRDKF-Open', ...
               'RDKFLOC', 'SRDKFLOC-Closed', 'SRDKFLOC-Open'};
  estFields = {'ckfSample', 'crkfSample', 'dseacpSample', 'dkfSample', 'rdkfSample', ...
               'sdkfClosedSample', 'sdkfOpenSample', 'srdkfClosedSample', 'srdkfOpenSample', ...
               'rdkflocSample', 'srdkflocClosedSample', 'srdkflocOpenSample'};
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
