function y = cumsum( x, dim )

%
% Basic argument check
%

s = x.size_;
switch nargin,
    case 0,
        error( 'Not enough input arguments.' );
    case 1,
        dim = cvx_default_dimension( s );
    case 2,
        if cvx_check_dimension( dim, false ),
            error( 'Second argument must be a dimension.' );
        end
end

if dim > length( s ) | s( dim ) <= 1,

    y = x;

else

    b = x.basis_;
    sb = size( b );
    need_perm = any( s( 1 : dim - 1 ) > 1 );
    if need_perm,
        ndxs = reshape( 1 : prod( s ), s );
        ndxs = permute( ndxs, [ dim, 1 : dim - 1, dim + 1 : length( s ) ] );
        b = b( :, ndxs );
    end
    b = reshape( b, prod( sb ) / s( dim ), s( dim ) );
    b = cumsum( b, 2 );
    b = reshape( b, sb );
    if need_perm,
        b( :, ndxs ) = b;
    end
    y = cvx( s, b );
    v = cvx_vexity( y );
    if any( isnan( v( : ) ) ),
        error( sprintf( 'Disciplined convex programming error:\n   Addition of convex and concave terms is forbidden.' ) );
    end

end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
