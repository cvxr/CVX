function y = transpose( x );
error( cvx_verify( x ) );

%
% Check sizes
%

s = x.size_;
if length( s ) > 2,
    error( 'Transpose of an ND array is not defined.' );
end

%
% Perform transpose
%

b = x.basis_;
sv = size( b );
[ r, c, v ] = find( b );
r = r - 1;
r = rem( r, s( 1 ) ) * s( 2 ) + floor( r / s( 1 ) ) + 1;
y = cvx( problem( x ), [ s( 2 ), s( 1 ) ], sparse( r, c, v, sv( 1 ), sv( 2 ) ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
