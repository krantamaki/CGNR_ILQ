% Updating implementation of trilsolve function used to solve
% linear system Lx = b, where L is invertible lower triangular matrix
function x = trilsolve(L, b)

  nrows = size(L, 1);
  x = b;

  for i=1:nrows

    x(i) = b(i) / L(i, i);

    b((i + 1):nrows) = b((i + 1):nrows) - L((i + 1):nrows, i) * x(i);

  end

end
