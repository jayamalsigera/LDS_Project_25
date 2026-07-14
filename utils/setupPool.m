function setupPool()
%SETUPPOOL  Size the parallel pool to the SLURM CPU allocation (if any).
%
%   Under SLURM the default 'Processes' profile sizes the pool to the node's
%   physical core count, which can exceed --cpus-per-task and oversubscribe
%   the allocation (and blow past the memory request, since peak memory scales
%   with the number of live workers). When SLURM_CPUS_PER_TASK is set this
%   starts a pool with exactly that many workers. Off the cluster (variable
%   unset) it is a no-op and parfor uses the default pool.

  n = str2double(getenv('SLURM_CPUS_PER_TASK'));
  if isnan(n) || n < 1
    return;   % not under SLURM -> leave default pool behaviour
  end
  p = gcp('nocreate');
  if isempty(p)
    parpool(n);
  elseif p.NumWorkers ~= n
    delete(p);
    parpool(n);
  end
end
