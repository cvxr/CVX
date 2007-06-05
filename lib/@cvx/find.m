function [ i, j, v ] = find( x )

ndxs = find( any( x.basis_, 1 ) );
ndxs = ndxs( : );

if nargout > 1,
    i = ndxs - 1;
    j = floor( i / x.size_(1) ) + 1;
    i = rem( i, x.size_(1) ) + 1;
end

if nargout > 2,
    v = reshape( cvx_subsref( x, ndxs ), length( i ) );
end

