function [plant, m] = sstPlant(P)
%SSTPLANT  Build the single-target tracking plant for a params struct.
%
%   [plant, m] = sstPlant(P) constructs the 2D or 3D plant from a params
%   struct P (as returned by collectParams), dispatching on P.dim, and returns
%   the per-sensor measurement dimension m (= P.dim). This is the one
%   dimension-specific step shared by every tune* function, so the sweep logic
%   itself stays identical across 2D and 3D.

  switch P.dim
    case 2
      plant = SingleTarget2dModel(P.Ts, P.sensorCount, P.outputNoiseStd, P.T, P.turnRate);
    case 3
      plant = SingleTarget3dModel(P.Ts, P.sensorCount, P.noiseScale, P.T);
    otherwise
      error('sstPlant:dim', 'Unsupported dim = %g (expected 2 or 3).', P.dim);
  end
  m = P.dim;   % measurement rows per sensor
end
