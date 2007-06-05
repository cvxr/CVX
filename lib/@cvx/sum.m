function y = sum( x, dim )

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
    y = cvx( s, sparse( 1, prod( s ) ) );

else

    p  = prod( s( 1 : dim - 1 ) );
    cc = 0 : prod( s ) - 1;
    cl = rem( cc, p );
    cr = floor( cc / ( p * s( dim ) ) );
    cc = cl + cr * p + 1;
    s( dim ) = 1;
    b = x.basis_;
    [ r, c, v ] = find( b );
    b = sparse( r, cc( c ), v, size( b, 1 ), prod( s ) );
    y = cvx( s, b );
    v = cvx_vexity( y );
    if any( isnan( v( : ) ) ),
        error( sprintf( 'Disciplined convex programming error:\n   Addition of convex and concave terms is forbidden.' ) );
    end

end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
