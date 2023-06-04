% Implementation of CGNR linear solver with ILQ preconditioning

function x = cgnr_ilq(A, x0, b, max_iter, tol, param)

    A_T = A';
    L = ilq(A_T, param);
    L = inv_L(L, param);
    % L = inv(L)
    L_T = L';
  
    x = x0;

    b_tmp = L * b;
    b_tmp = L_T * b_tmp;
    b_tmp = A_T * b_tmp;

    r = b_tmp - ilq_mult(A_T, L_T, L, A, x);
    p = r;

    for iter = 1:max_iter

        p_tmp = ilq_mult(A_T, L_T, L, A, p);
      
        alpha = dot(p, r) / dot(p, p_tmp);

	x = x + alpha * p;

	r_tmp = r;

	r = r - alpha * p_tmp;
	p = r + dot(r, r) / dot(r_tmp, r_tmp) * p;

	if (norm(r) < tol)
	    disp(["Required number of iterations: ", num2str(iter)])
            return;
	end
    end
    
    disp("WARNING: Did not converge to wanted tolerance within given number of iterations!")
    
end
