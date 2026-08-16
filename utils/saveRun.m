%% Persist a Monte Carlo run to results/<file>.mat.
%
% Builds a filename that encodes the most commonly varied parameters
% (T, nodeCount, sensorCount, klTolerance, totalRuns) plus a timestamp,
% and stores the full parameter set inside the .mat alongside the results.
%
function savedPath = saveRun(scriptName, params, extras, netGraph, results, samples)

  ts = char(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
  totalRuns = extras.totalRuns;

  isLfmScript = contains(scriptName, 'Lfm');
  if isLfmScript && isfield(params, 'lfmKlTolerance')
    bTag = sprintf('b%g_blfm%g', params.klTolerance, params.lfmKlTolerance);
  else
    bTag = sprintf('b%g', params.klTolerance);
  end

  tag = sprintf('T%d_N%ds%d_%s_runs%d', ...
                params.T, params.nodeCount, params.sensorCount, ...
                bTag, totalRuns);

  outDir = fullfile(projectRoot(), 'results');
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
