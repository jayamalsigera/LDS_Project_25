function plotSST3dTrajectories(path)
%PLOTSST3DTRAJECTORIES  Plot truth and estimated 3D trajectories for a saved run.
%
% Position lives in state rows 4:6 ([px py pz]), so trajectories are drawn with
% plot3. Works for both the
% nominal-plant script (estimateSST3d) and the least-favorable model script
% (estimateSST3dLfm); the LFM branch is taken when samples has a nomSample field.

  runData = loadRun(path);

  samples = runData.samples;
  isLfm   = isfield(samples, 'nomSample');

  % Truth trajectory
  if isLfm
    figure
    plot3(samples.nomSample.X(4, :), samples.nomSample.X(5, :), ...
          samples.nomSample.X(6, :), 'DisplayName', 'Nominal');
    hold on
    plot3(samples.mdlSample.X(4, :), samples.mdlSample.X(5, :), ...
          samples.mdlSample.X(6, :), 'DisplayName', 'LFM');
    hold off
    title("Simulated Trajectory")
    xlabel('$p_x$', 'Interpreter', 'latex');
    ylabel('$p_y$', 'Interpreter', 'latex');
    zlabel('$p_z$', 'Interpreter', 'latex');
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
      plotEstVsTruth(est, estLabels{k}, samples.mdlSample.X);
    end
  end
end

function meanX_hat = estMean(estSample)
  if isprop(estSample, 'X_hat')
    meanX_hat = squeeze(mean(estSample.X_hat, 2));
  else
    meanX_hat = estSample.x_hat;
  end
end

function plotEstVsTruth(estSample, label, Xtruth)
  meanX_hat = estMean(estSample);
  figure
  plot3(meanX_hat(4, :), meanX_hat(5, :), meanX_hat(6, :), 'DisplayName', label);
  hold on
  plot3(Xtruth(4, :), Xtruth(5, :), Xtruth(6, :), 'DisplayName', 'Truth');
  hold off
  title(sprintf("%s Estimated Trajectory", label));
  xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
  ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
  zlabel('$\hat{p}_z$', 'Interpreter', 'latex');
  legend(); grid();
end

function plotEstVsTruthWithNominal(estSample, label, Xlfm, Xnom)
  meanX_hat = estMean(estSample);
  figure
  plot3(meanX_hat(4, :), meanX_hat(5, :), meanX_hat(6, :), 'DisplayName', label);
  hold on
  plot3(Xlfm(4, :), Xlfm(5, :), Xlfm(6, :), 'DisplayName', 'LFM truth');
  plot3(Xnom(4, :), Xnom(5, :), Xnom(6, :), '--', 'DisplayName', 'Nominal');
  hold off
  title(sprintf("%s Estimated Trajectory", label));
  xlabel('$\hat{p}_x$', 'Interpreter', 'latex');
  ylabel('$\hat{p}_y$', 'Interpreter', 'latex');
  zlabel('$\hat{p}_z$', 'Interpreter', 'latex');
  legend(); grid();
end
