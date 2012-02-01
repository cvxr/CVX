function cvx_optval = sigma_max( x )

%SIGMA_MAX   Internal cvx version.

error( nargchk( 1, 1, nargin ) );
if ndims( x ) > 2,
    error( 'lambda_max is not defined for N-D arrays.' );
elseif ~cvx_isaffine( x ),
    error( 'Input must be affine.' );
end

%
% Construct problem
% 

[ m, n ] = size( x );
cvx_optval = lambda_max( [ zeros( m, m ), x ; x', zeros( n, n ) ] );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
