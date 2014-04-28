function bx = cvx_sparsify( bx, tt, mode )

global cvx___
if isempty( tt ),
    by = bx;
elseif ~any( tt ),
    return
else
    by = bx( :, tt );
end
[ xR, by ] = cvx_bcompress( by, mode, 0 );
nx = size( by, 2 );
ndim = cvx_pushvar( nx, cvx_classify_mex( by, cvx___.classes ) );
bz = sparse( ndim, 1 : nx, 1 );
nz = size(bz,1);
by( nz, end ) = 0;
cvx_pushcnstr( by - bz, true );
bz = bz * xR;
if isempty( tt ),
    bx = bz;
else
    bx( 1:nz, tt ) = bz;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
