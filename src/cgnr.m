% Implemetation of the conjugate gradient on the normal equations method for
% benchmarking purposes

function x = cgnr(A, x0, b, max_iter, tol)

  x = x0;
  b_tmp = (b' * A)';

  x_tmp = A * x0;
  x_tmp = (x_tmp' * A)';
  
  r = b_tmp - x_tmp;
  p = r;

  for iter = 1:max_iter

    disp(["Iteration: ", num2str(iter)])

    p_tmp = A * p;
    p_tmp = (p_tmp' * A)';

    alpha = (p' * r) / (p' * p_tmp);
    x = x + alpha * p;

    r_tmp = r;

    r = r - alpha * p_tmp;
    p = r + (r' * r) / (r_tmp' * r_tmp) * p;

    if (norm(r) < tol)
      disp(["Required number of iterations: ", num2str(iter)])
      return;
    end
  end

  disp("WARNING: Did not converge to wanted tolerance within given number of iterations!")
  
end
