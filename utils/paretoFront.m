function isFront = paretoFront(x, y)
%PARETOFRONT  Non-dominated (minimising) subset of a 2D point cloud.
%
%   isFront = paretoFront(x, y) returns a logical column vector marking the
%   points (x_i, y_i) that are Pareto-optimal when both axes are minimised
%   (here: transmission rate and RMSE). Point i is dominated when some other
%   point j is no worse on both axes and strictly better on at least one.

  x = x(:);
  y = y(:);
  n = numel(x);
  isFront = true(n, 1);
  for i = 1:n
    dominated = (x <= x(i)) & (y <= y(i)) & ((x < x(i)) | (y < y(i)));
    dominated(i) = false;
    if any(dominated)
      isFront(i) = false;
    end
  end
end
