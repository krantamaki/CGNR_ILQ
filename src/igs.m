% Implementation of incomplete Gram-Schmidt process
% i.e. only does the projections to the p preceeding vectors
% NOTE! Does the process for the transpose of A i.e. orthogonalises
% the rows rather than the columns

function E = igs(A, p)

    E = A;
    nrows = size(A, 1);

    for k=1:p:nrows

      p0 = min([nrows - k, p]);

      E(k, :) = E(k, :) / norm(E(k, :));
    
      for i = 1:p0 - 1

	for j = 1:i

	  E(k + i, :) = E(k + i, :) - proj(A(k + i, :), E(k + i - j, :));
	  
	end

	E(k, :) = E(k, :) / norm(E(k, :));
	
      end

    end
    
end

