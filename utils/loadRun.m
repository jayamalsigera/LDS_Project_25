%% Load a saved Monte Carlo run produced by saveRun.
function runData = loadRun(path)

  s = load(path);
  if ~isfield(s, 'runData')
    error('loadRun:badFile', ...
          '%s does not contain a ''runData'' variable', path);
  end
  runData = s.runData;
end
