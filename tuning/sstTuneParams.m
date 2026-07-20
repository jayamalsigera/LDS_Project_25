function P = sstTuneParams(scenario)
%SSTTUNEPARAMS  Validate a tune* scenario token and load its params struct.
%
%   P = sstTuneParams(scenario) checks that SCENARIO is 'sst2d' or 'sst3d'
%   (required -- there is no default, so a 2D/3D mix-up fails loudly) and
%   returns the corresponding params struct via collectParams. This is the
%   shared front door for every tune* function.

  valid = {'sst2d', 'sst3d'};
  if nargin < 1 || ~(ischar(scenario) || isstring(scenario)) ...
      || ~ismember(char(scenario), valid)
    error('sstTuneParams:scenario', ...
          "Pass a scenario: 'sst2d' or 'sst3d' (required, no default).");
  end
  P = collectParams([char(scenario) 'Params']);
end
