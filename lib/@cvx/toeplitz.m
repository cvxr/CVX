function t = toeplitz( c, r )
error( nargchk( 1, 2, nargin ) );

%
% Check arguments
%

error( cvx_verify( c ) );
if nargin < 2,
    prob = problem( c );
    c = cvx_subasgn( c, 1, cvx_subref( c, 1 ) );
    r = c;
    c = conj( c );
else,
    error( cvx_verify( r ) );
    [ prob, c, r ] = cvx_operate( [], c, r );
    temp = cvx_subref( r, 1 ) - cvx_subref( c, 1 );
    if ~cvx_isconstant( temp ) | cvx_constant( temp ) ~= 0,
        warning('MATLAB:toeplitz:DiagonalConflict',['First element of ' ...
               'input column does not match first element of input row. ' ...
               '\n         Column wins diagonal conflict.'])
    end
end

%
% Compute indices and construct data vector
% 

r = vec( r );
c = vec( c );
p = length( r );
m = length( c );
x = [ cvx_subref( r, p : -1 : 2 ) ; c ];

%
% Construct matrix
%

cidx = [ 0 : m - 1 ]';
ridx = p : -1 : 1;
t    = cidx( :, ones( p, 1 ) ) + ridx( ones( m, 1 ) , : );
t    = reshape( cvx_subref( x, t( : ) ), size( t ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
