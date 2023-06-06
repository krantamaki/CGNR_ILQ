% Computes the incomplete LQ decomposition of a given matrix
% with independent blocks of rows

function [L, Q] = bmilq(A, p)

   [m, n] = size(A);

   Q = sparse(zeros(m, n));
   L = sparse(zeros(n, n));

   for k = 1:p:m

     p0 = min([m - k, p]);

     for i = k:k + p0

       q = A(i, :);

       for j = k:i - 1

	 L(i, j) = q * Q(j, :)';
	 q = q - L(i, j) * Q(j, :);

       end

       L(i, i) = norm(q);
       Q(i, :) = q / L(i, i);

   end

end
