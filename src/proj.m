% Projection operator associated with Gram-Schmidt process

function p = proj(a, u)

    p = (dot(a, u) / dot(u, u)) * u;

end
