% Implementation of CGNR linear solver with ILQ preconditioning

function x = cgnr_ilq(A, x0, b, max_iter, tol, param, comp_cond)
  
    [L, Q] = milq(A, param);

    if comp_cond
      disp("\nComputing the condition numbers. This will take a while...\n")

      condA = cond(A);

      disp(["Condition number of A: ", num2str(condA)])

      disp(["Condition number of system of normal equations: ", num2str(condA ^ 2)])

      L_inv = inv(L);
      
      disp(["Condition number of preconditioned system: ", num2str(cond(A' * L_inv' * L_inv * A)), "\n"])

    end
  
    x = x0;

    b_tmp = trilsolve(L, b);
    b_tmp = triulsolve(L, b_tmp);
    b_tmp = (b_tmp' * A)';

    r = b_tmp - ilq_mult(A, L, x);
    p = r;

    for iter = 1:max_iter

        disp(["Iteration: ", num2str(iter)])

        p_tmp = ilq_mult(A, L, p);
      
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
