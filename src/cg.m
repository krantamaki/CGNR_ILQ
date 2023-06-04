% Classic implemetation of the conjugate gradient method for
% benchmarking purposes

function x = cg(A, x0, b, max_iter, tol)

  x = x0;
  r = b - A * x0;
  p = r;

  for iter = 1:max_iter

    alpha = (p' * r) / (p' * A * p);
    x = x + alpha * p;

    r_tmp = r;

    r = r - alpha * A * p;
    p = r + (r' * r) / (r_tmp' * r_tmp) * p;

    if (norm(r) < tol)
      disp(["Required number of iterations: ", num2str(iter)])
      return;
    end
  end

  disp("WARNING: Did not converge to wanted tolerance within given number of iterations!")
  
end
