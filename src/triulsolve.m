% Implementation of triusolve function used to solve
% linear system Ux = b, but where U is given as L^T.
% Thus, to avoid transposing a matrix this function solves
% related problem x^T L = b^T
function x = triulsolve(L, b)

  nrows = size(L, 1);
  b_T = b';
  x_T = b_T;

  for i=nrows:-1:1

    x_T(i) = b_T(i) / L(i, i);

    b_T(1:(i-1)) = b_T(1:(i-1)) - L(i, 1:(i-1)) * x_T(i);

  end

  x = x_T';

end
