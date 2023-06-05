% Implementation of the incomplete LQ decomposition for given matrix A

function L = ilq(A, p)

    nrows = size(A, 1);
    L = sparse(zeros(size(A)));
    
    disp("Doing the Incomplete Gram-Schmidt process")
    
    Q = igs(A, p);

    disp("Doing the Incomplete LQ decomposition")

    L(1, 1) = dot(Q(:, 1), A(1, :));
    
    for i = 2:nrows

        p0 = min([i, p]);

	L(i, i) = dot(Q(:, i), A(i, :));

	for j = 1:p0 - 1

	    L(i, i - j) = dot(Q(:, j), A(i, :));
	  
	end
    end
end
