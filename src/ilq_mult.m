% Function that does the required A_T * L_T * L * A * x multiplication
% as a chain of matrix vector multiplications rather than matrix matrix
% multiplications

function x = ilq_mult(A_T, L_T, L, A, x0)
 
    x = x0;

    x = A * x;
    x = L * x;
    x = L_T * x;
    x = A_T * x;

end
