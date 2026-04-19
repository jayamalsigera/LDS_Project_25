function p = collectParams()
%COLLECTPARAMS  Return a struct of all variables declared in sst2dParams.m.
%
% Keeps sst2dParams.m as the single source of truth for shared simulation
% parameters; this helper packages them into a struct for logging/saving.

  sst2dParams;
  names = who;
  p = struct();
  for k = 1:numel(names)
    if strcmp(names{k}, 'names')
      continue;
    end
    p.(names{k}) = eval(names{k});
  end
end
