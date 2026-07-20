function p = collectParams(paramsScript)
%COLLECTPARAMS  Return a struct of all variables declared in a params script.
%
% Keeps the shared params files (sst2dParams.m / sst3dParams.m) as the single
% source of truth for simulation parameters; this helper sources one of them and
% packages the resulting variables into a struct for logging/saving.
%
%   collectParams()               -> sources sst2dParams (default, 2D)
%   collectParams('sst2dParams')  -> 2D
%   collectParams('sst3dParams')  -> 3D

  if nargin < 1 || isempty(paramsScript)
    paramsScript = 'sst2dParams';
  end

  eval(paramsScript);           % populate this workspace from the params file
  names = who;
  p = struct();
  for k = 1:numel(names)
    if any(strcmp(names{k}, {'names', 'paramsScript'}))
      continue;
    end
    p.(names{k}) = eval(names{k});
  end
end
