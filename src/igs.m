% Implementation of incomplete Gram-Schmidt process
% i.e. only does the projections to the p preceeding vectors

function E = igs(A, p)

    E = zeros(size(A));
    ncols = size(A, 2);

    E(:, 1) = A(:, 1) / norm(A(:, 1));

    for k=2:ncols

        p0 = min([k, p]);
        E(:, k) = A(:, k);
    
        for i = 1:p0 - 1
            E(:, k) = E(:, k) - proj(A(:, k), E(:, k - i));
        end

        E(:, k) = E(:, k) / norm(E(:, k));

    end
end
