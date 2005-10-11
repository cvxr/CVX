function z = kron( x, y )
error( nargchk( 2, 2, nargin ) );
error( cvx_verify( x, y ) );
[ prob, x, y ] = cvx_operate( [], x, y );

%
% Enforce DCP rules
%

xc = cvx_isconstant( x );
yc = cvx_isconstant( y );
if ~xc & ~yc,
    error( 'The multiplication of non-cvx_constant expressions is forbidden.' );
end

%
% Determine new size
%

sx = size( x );
sy = size( y );
dx = length( sx );
dy = length( sy );
dz = max( dx, dy );
sx = [ sx, ones( 1, dz - dx ) ];
sy = [ sy, ones( 1, dz - dy ) ];
sz = sx .* sy;

%
% Construct indices
%

for k = 1 : dz,
    temp = 0 : sz( k ) - 1;
    ndxx{k} = floor( temp / sy( k ) ) + 1;
    ndxy{k} = rem( temp, sy( k ) ) + 1;
end
temp = reshape( 1 : prod( sx ), sx );
ndxx = temp( ndxx{:} );
temp = reshape( 1 : prod( sy ), sy );
ndxy = temp( ndxy{:} );

%
% Expand and multiply
%

bx = cvx_basis( x );
by = cvx_basis( y );
if xc, 
    bz = diag( bx( ndxx, 1 ) ) * by( ndxy, : );
else,
    bz = diag( by( ndxy, 1 ) ) * bx( ndxx, : );
end

%
% Create object
%

z = cvx( prob, sz, bz );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
