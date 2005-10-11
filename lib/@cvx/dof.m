function y = dof( cv )
error( cvx_verify( cv ) );

x = cvx_basis( cv );
if ~isreal( x ),
    x = [ real( x ) ; imag( x ) ];
end

%
% Extract the non-zero rows and columns
%

x = x( any( x, 2 ), any( x, 1 ) );
[ m, n ] = size( x );

%
% Normalize rows by the sign of the first non-zero column
%

[ c, r, v ] = find( x' );
c = c( : ); r = r( : ); v = v( : );
temp = r( v < 0 & [ 1 ; diff( r ) ] ~= 0 );
x( temp, : ) = - x( temp, : );

%
% Sort the rows and count the unique ones
%

ndxs = 1 : m;
for k = n : -1 : 1,
   [ v, ind ] = sort( x( ndxs, k ) );
   ndxs = ndxs( ind );
end
x = x( ndxs, : );
d = x( 2 : end, : ) ~= x( 1 : end - 1, : );
d = [ 1 ~= 0; any( d, 2 ) ];

y = full( sum( d ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
