function savedPath = saveRun(scriptName, params, extras, netGraph, results, samples)
%SAVERUN  Persist a Monte Carlo run to results/<scriptName>/<file>.mat.
%
% Builds a filename that encodes the most commonly varied parameters
% (T, nodeCount, sensorCount, klTolerance, totalRuns) plus a timestamp,
% and stores the full parameter set inside the .mat alongside the results.

  ts = char(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
  totalRuns = extras.totalRuns;
  tag = sprintf('T%d_N%ds%d_b%g_runs%d', ...
                params.T, params.nodeCount, params.sensorCount, ...
                params.klTolerance, totalRuns);

  outDir = fullfile(projectRoot(), 'results', scriptName);
  if ~exist(outDir, 'dir')
    mkdir(outDir);
  end

  fname = sprintf('%s_%s_%s.mat', scriptName, tag, ts);
  savedPath = fullfile(outDir, fname);

  runData = struct();
  runData.script = scriptName;
  runData.timestamp = ts;
  runData.gitSha = tryGitSha(projectRoot());
  runData.params = params;
  runData.extras = extras;
  runData.netGraph = netGraph;
  runData.results = results;
  runData.samples = samples;

  save(savedPath, 'runData', '-v7.3');
  fprintf('Saved run to %s\n', savedPath);
end

function sha = tryGitSha(repoRoot)
  cmd = sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null', repoRoot);
  [status, out] = system(cmd);
  if status == 0
    sha = strtrim(out);
  else
    sha = '';
  end
end

function root = projectRoot()
  root = fileparts(fileparts(mfilename('fullpath')));
end
