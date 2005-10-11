function y = blkdiag( varargin )
error( cvx_verify( varargin{:} ) );
[ prob, varargin{ : } ] = cvx_operate( [], varargin{ : } );

nv = 0;
nz = 0;
sz = [ 0, 0 ];
for k = 1 : nargin,
    x  = varargin{k};
    sx = size( x );
    if length( sx ) > 2,
        error( 'N-D matrices not supported.' );
    end
    b  = cvx_basis( x );
    sz = sz + sx;
    nv = max( nv, size( b, 2 ) );
    nz = nz + nnz( b );
end
bz = sparse( [], [], [], prod( sz ), nv, nz );
roff = 0;
coff = 0;
for k = 1 : nargin,
    x  = varargin{k};
    b  = cvx_basis( x );
    sx = size( x );
    ndxr = [ roff : roff + sx( 1 ) - 1 ]';
    ndxr = ndxr( :, ones( 1, sx( 2 ) ) );
    ndxc = [ coff : coff + sx( 2 ) - 1 ];
    ndxc = ndxc( ones( 1, sx( 1 ) ), : );
    bz( ndxc( : ) * sz( 1 ) + ndxr( : ) + 1, 1 : size( b, 2 ) ) = b;
    roff = roff + sx( 1 );
    coff = coff + sx( 2 );
end
y = cvx( prob, sz, bz );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
