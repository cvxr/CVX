function y = diag( v, k )
error( cvx_verify( v ) );

switch nargin,
    case 0,
        error( 'Not enough arguments.' );
    case 1,
        k = 0;
    case 2,
        if ~isnumeric( k ) | k ~= floor( k ),
            error( 'Second argument must be an integer.' );
        end
end
 
s = size( v );

if length( s ) ~= 2,
    error( 'First input must be 2D.' );
end

if k < 0,
    absk = -k;
    roff = absk;
    coff = 0;
else,
    absk = +k;
    roff = 0;
    coff = absk;
end

if any( s == 1 ),
    [ r, c, d ] = find( cvx_basis( v ) );
    r = r + roff + ( r + coff - 1 ) * ( s( 1 ) + absk ); 
    s = s( [ 1, 1 ] ) + absk;
    y = cvx( problem( v ), s, sparse( r, c, d, prod( s ), max( c ) ) );
elseif -k >= s( 1 ) | k >= s( 2 ),
    y = cvx( problem( v ), [ 0, 0 ], sparse( 0, 0 ) );
else,
    [ r, c, d ] = find( cvx_basis( v ) );
    r = r - 1;
    rows = rem( r, s( 1 ) );
    cols = floor( r / s( 1 ) );
    temp = rows - roff == cols - coff;
    r = rows( temp ) - roff + 1;
    c = c( temp );
    d = d( temp );
    y = cvx( problem( v ), [ length( r ), 1 ], sparse( r, c, d, length( r ), max( c ) ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
