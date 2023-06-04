% Efficient implementation of Gaussian elimination for inverting
% lower diagonal matrix L with p subdiagonals filled

function L_inv = inv_L(L, p)

    nrows = size(L, 1);
    L_inv = eye(size(L));
    
    for i = 1:nrows

        L_inv(i, :) = L_inv(i, :) / L(i, i);

        p0 = min([p, nrows - i]);

	for j = 1:p0
	  
	    L_inv(i + j, :) = L_inv(i + j, :) - L_inv(i, :) * L(i + j, i);
	  
	end
    end
end
