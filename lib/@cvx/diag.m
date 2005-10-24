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
    y = cvx( problem( v ), [ 0, 1 ], sparse( 0, 0 ) );
else,
    b = cvx_basis( v );
    sz = min( s - [ roff, coff ] );
    r = ( roff + s( 1 ) * coff + 1 ) + ( 0 : sz - 1 ) * ( s( 1 ) + 1 );
    y = cvx( problem( v ), length( r ),  b( r, : ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
