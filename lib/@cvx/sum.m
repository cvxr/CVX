function y = sum( x, dim )
error( cvx_verify( x ) );

%
% Basic argument check
%

s = size( x );
switch nargin,
    case 0,
        error( 'Not enough input arguments.' );
    case 1,
        dim = [ find( s > 1 ), 1 ];
        dim = dim( 1 );
    case 2,
        if ~isnumeric( dim ) | dim <= 0 | dim ~= floor( dim ),
            error( 'Second argument must be a dimension.' );
        end
end

if dim > length( s ) | s( dim ) == 1,

    y = x;
    
elseif s( dim ) == 0,

    s( dim ) = 1;
    y = cvx( problem( x ), s, sparse( prod( s ), 0 ) );
    
else,

    b = cvx_basis( x );
    [ r, c, v ] = find( b );
    p  = prod( s( 1 : dim - 1 ) );
    r  = r - 1;
    rl = rem( r, p );
    rr = floor( r / ( p * s( dim ) ) );
    r  = rl + rr * p + 1;
    s( dim ) = 1;
    y = cvx( problem( x ), s, sparse( r, c, v, prod( s ), size( b, 2 ) ) );
    v = cvx_vexity( y );
    if any( isnan( v( : ) ) ),
        error( sprintf( 'Disciplined convex programming error:\n   Addition of convex and concave terms is forbidden.' ) );
    end

end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
