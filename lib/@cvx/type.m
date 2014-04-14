function st = type( x )

s = size( x );
len = prod( s );
if len == 1,
    nzs = len;
    st = 'scalar';
else
    nzs = nnz( any( cvx_basis( x ), 1 ) );
    nd = length( s );
    st = sprintf( '%dx', s );
    st = st( 1 : end - 1 );
    if nd > 2,
        st = [ st, ' array' ];
    elseif any( s == len ),
        st = [ st, ' vector' ];
    else
        st = [ st, ' matrix' ];
    end
end
if nzs < len,
    if nzs > 1,
        st = sprintf( '%s, %d nonzeros', st, nzs );
    elseif nzs == 1,
        st = sprintf( '%s, 1 nonzero', st );
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
