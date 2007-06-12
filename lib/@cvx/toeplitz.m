function t = toeplitz( c, r )
error( nargchk( 1, 2, nargin ) );

%
% Check arguments
%

if nargin < 2,
    c    = vec( c );
    m    = length( c );
    p    = m;
    c    = [ cvx_subsref( c, p : -1 : 1 ) ; conj( cvx_subsref( c, 2 : p ) ) ]
else
    temp = cvx_subsref( r, 1 ) - cvx_subsref( c, 1 );
    if ~cvx_isconstant( temp ) | cvx_constant( temp ) ~= 0,
        warning('MATLAB:toeplitz:DiagonalConflict',['First element of ' ...
               'input column does not match first element of input row. ' ...
               '\n         Column wins diagonal conflict.'])
    end
    r = vec( r );
    c = vec( c );
    p = length( r );
    m = length( c );
    x = [ cvx_subsref( r, p : -1 : 2 ) ; c ];
end

%
% Construct matrix
%

cidx = [ 0 : m - 1 ]';
ridx = p : -1 : 1;
t    = cidx( :, ones( p, 1 ) ) + ridx( ones( m, 1 ) , : );
t    = reshape( cvx_subsref( x, t ), size( t ) );

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
