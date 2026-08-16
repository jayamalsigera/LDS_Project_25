%% Return a struct of all variables declared in sst3dParams.m.
%
% Keeps sst3dParams.m as the single source of truth for simulation parameters;
% this helper sources it and packages the resulting variables into a struct for
% logging/saving.
%
function p = collectParams()

  sst3dParams;                  % populate this workspace from the params file
  names = who;
  p = struct();
  for k = 1:numel(names)
    if strcmp(names{k}, 'names')
      continue;
    end
    p.(names{k}) = eval(names{k});
  end
end
