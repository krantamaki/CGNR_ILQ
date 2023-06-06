% Computes the incomplete LQ decomposition of a given matrix

function [L, Q] = milq(A, p)

   [m, n] = size(A);

   Q = zeros(m, n);
   L = zeros(n, n);

   for i = 1:m

     p0 = min([m - i, p]);

     q = A(i, :);

     for j = i:i + p0

       L(i, j) = q * Q(j, :)';
       q = q - L(i, j) * Q(j, :);
      
     end

     L(i, i) = norm(q);
     Q(i, :) = q / L(i, i);

   end

end
