% Projection operator associated with Gram-Schmidt process

function p = proj(a, u)

    p = ((a * u') / (u * u')) * u;

end
