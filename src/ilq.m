% Implementation of the incomplete LQ decomposition for given matrix A

function L = ilq(A_T, p)

    ncols = size(A_T, 2);
    L = zeros(size(A_T));
  
    Q = igs(A_T, p);

    L(1, 1) = dot(Q(:, 1), A_T(:, 1));
    
    for i = 2:ncols

        p0 = min([i, p]);

	L(i, i) = dot(Q(:, i), A_T(:, i));

	for j = 1:p0 - 1

	    L(i, i - j) = dot(Q(:, j), A_T(:, i));
	  
	end
    end
end
