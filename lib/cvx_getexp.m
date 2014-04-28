function by = cvx_getexp( bx )

global cvx___
nb = size( bx, 2 );
if nb == 0,
    by = bx;
    return; 
end
cc = bx( 1, : );
if nnz( cc ),
    bx( 1, : ) = 0;
    cc = exp( cc );
else
    cc = 1;
end
ry = ones( nb, 1 );
if nnz( bx ),
    bz = sum( bx ~= 0, 1 );
    tt = ( bz > 1 ) | ( ( bz > 0 ) & ( sum( bx, 1 ) ~= 1 ) );
    bx = cvx_sparsify( bx, tt, 'none' );
    [ rx, cx ] = find( bx );
    if ~isempty( rx ),
        cvx___.exponential( end + 1 : max(rx), 1 ) = 0;
        exps = cvx___.exponential( rx, 1 );
        tt = exps == 0;
        if any( tt ),
            exps(tt) = cvx_pushexp( rx(tt) );
        end
        ry(cx) = exps;
    end
end
by = sparse( ry, 1:nb, cc );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
