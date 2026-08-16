%% Load a saved Monte Carlo run produced by saveRun.
%
% With no argument, loads the most recently modified .mat in results/.
%
function runData = loadRun(path)
  resultsDir = fullfile(projectRoot(), 'results');

  if nargin < 1 || isempty(path)
    path = latestRun(resultsDir);
  else
    path = fullfile(resultsDir, path);
  end

  s = load(path);
  if ~isfield(s, 'runData')
    error('loadRun:badFile', ...
          '%s does not contain a ''runData'' variable', path);
  end
  runData = s.runData;
end

function path = latestRun(resultsDir)
  files = dir(fullfile(resultsDir, '*.mat'));
  if isempty(files)
    error('loadRun:noRuns', 'No .mat runs found in %s', resultsDir);
  end
  [~, idx] = max([files.datenum]);
  path = fullfile(resultsDir, files(idx).name);
  fprintf('Loading latest run %s\n', files(idx).name);
end

function root = projectRoot()
  root = fileparts(fileparts(mfilename('fullpath')));
end
