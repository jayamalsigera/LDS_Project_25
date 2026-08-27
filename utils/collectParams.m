%% Return a struct of all variables declared in sst3dParams.m.
%
% Keeps sst3dParams.m as the single source of truth for simulation parameters;
% this helper sources it and packages the resulting variables into a struct for
% logging/saving.
%
function p = collectParams(varargin)

  sst3dParams;                  % populate this workspace from the params file
  names = who;
  p = struct();
  for k = 1:numel(names)
    if strcmp(names{k}, 'names') || strcmp(names{k}, 'varargin')
      continue;
    end
    p.(names{k}) = eval(names{k});
  end

  if nargin == 0
    return;
  end

  if numel(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
  end

  i = 1;
  while i <= numel(varargin)
    arg = varargin{i};
    if isstruct(arg)
      f = fieldnames(arg);
      for j = 1:numel(f)
        p.(f{j}) = arg.(f{j});
      end
      i = i + 1;
    elseif (ischar(arg) || isstring(arg)) && contains(arg, '=')
      parts = split(string(arg), '=');
      k = char(strtrim(parts(1)));
      vStr = char(strtrim(join(parts(2:end), '=')));
      val = str2num(vStr); %#ok<ST2NM>
      if isempty(val)
        if strcmpi(vStr, 'true')
          val = true;
        elseif strcmpi(vStr, 'false')
          val = false;
        else
          val = vStr;
        end
      end
      p.(k) = val;
      i = i + 1;
    elseif i + 1 <= numel(varargin) && (ischar(arg) || isstring(arg))
      k = char(arg);
      val = varargin{i + 1};
      if (ischar(val) || isstring(val))
        numVal = str2num(char(val)); %#ok<ST2NM>
        if ~isempty(numVal)
          val = numVal;
        elseif strcmpi(char(val), 'true')
          val = true;
        elseif strcmpi(char(val), 'false')
          val = false;
        else
          val = char(val);
        end
      end
      p.(k) = val;
      i = i + 2;
    else
      i = i + 1;
    end
  end
end
