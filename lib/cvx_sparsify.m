function [ bx, ndim ] = cvx_sparsify( bx, tt, mode )

global cvx___
ox = isa( bx, 'cvx' );
if ox,
    sx = size( bx );
    bx = cvx_basis( bx );
end
if isempty( tt ),
    by = bx;
elseif ~any( tt ),
    return
else
    by = bx( :, tt );
end
if isempty( mode ),
    by = bx;
else
    [ xR, by ] = cvx_bcompress( by, mode, 0 );
end
nx = size( by, 2 );
ndim = cvx_newvar( nx, cvx_classify_mex( by, cvx___.classes ) );
bz = sparse( ndim, 1 : nx, 1 );
nz = size(bz,1);
by( nz, end ) = 0;
cvx_newcnstr( by - bz, true );
if ~isempty( mode ),
    bz = bz * xR;
end
if isempty( tt ),
    bx = bz;
else
    bx( 1:nz, tt ) = bz;
end
if ox,
    bx = cvx( sx, bx );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
