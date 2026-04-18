%% Root-finding by bisection.
%
% Solves fun(x) = 0 on the interval [a, b] by bisection, assuming fun is
% monotonically increasing with fun(a) <= 0 <= fun(b).
%
% Returns when either |fun(x)| < tol or `max_iter` iterations have elapsed.
%
function x = bisect(fun, a, b, tol, max_iter)
  for k = 1:max_iter
    x = 0.5 * (a + b);
    fx = fun(x);

    if abs(fx) < tol
      return
    end

    if fx < 0
      a = x;
    else
      b = x;
    end
  end
end
