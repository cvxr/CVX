function y = cumsum( x, dim )
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

if dim > length( s ) | s( dim ) <= 1,

    y = x;
    
else,
    
    b = cvx_basis( x );
    sb = size( b );
    need_perm = any( s( 1 : dim - 1 ) > 1 );
    if need_perm,
        ndxs = reshape( 1 : prod( s ), s );
        ndxs = permute( ndxs, [ dim, 1 : dim - 1, dim + 1 : length( s ) ] );
        b = b( ndxs, : );
    end
    b = reshape( b, s( dim ), prod( sb ) / s( dim ) );
    b = cumsum( b, 1 );
    b = reshape( b, sb );
    if need_perm,
        b( ndxs, : ) = b;
    end
    y = cvx( problem( x ), s, b );
    v = cvx_vexity( y );
    if any( isnan( v( : ) ) ),
        error( sprintf( 'Disciplined convex programming error:\n   Addition of convex and concave terms is forbidden.' ) );
    end

end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
