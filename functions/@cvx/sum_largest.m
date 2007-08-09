function cvx_optval = sum_largest( x, k, dim )

%SUM_LARGEST   Internal cvx version.

error( nargchk( 2, 3, nargin ) );
sx = size( x );
if nargin < 3 | isempty( dim ),
    dim = [ find( sx > 1 ), 1 ];
    dim = dim( 1 );
elseif ~isnumeric( dim ) | ~isreal( dim ) | dim <= 0 | dim ~= floor( dim ),
    error( 'Second argument must be a dimension.' );
elseif ~isnumeric( k ) | ~isreal( k ),
    error( 'Third argument must be a real scalar.' );
elseif ~isreal( x ),
    error( 'First argument must be real.' );
end

%
% Quick exit cases
%

nd = max( dim, length( sx ) );
sx = [ sx, ones( 1, dim - nd ) ];
sy = sx;
sy( dim ) = 1;
if k <= 0,

    cvx_optval = zeros( sy );

elseif k <= 1,

    cvx_optval = k * max( x, [], dim );

elseif k >= sx( dim ),

    cvx_optval = sum( x, dim );

else

    nd = length( sx );
    ndxs = cell( 1, nd );
    [ ndxs{:} ] = deal( ':' );
    ndxs{ dim } = ones( 1, sx( dim ) );

    cvx_begin
        epigraph variable z( sy )
        variables xp( sx ) yp( sy )
        z == sum( xp, dim ) - k * yp;
        xp >= cvx_subsref( yp, ndxs{:} ) + x;
        xp >= 0;
    cvx_end

end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
