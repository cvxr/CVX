function cvx_optval = sigma_max( x )
error( nargchk( 1, 1, nargin ) );

%SIGMA_MAX    Maximum singular value.
%
%   SIGMA_MAX(X) returns the maximum singular value of X. X must be a 2-D
%   matrix, real or complex. SIGMA_MAX(X) is synonymous with NORM(X).
%
%   Disciplined quadratic programming information:
%       SIGMA_MAX(X) is convex and nonmontonic in X, so X must be affine.

%
% Quick exit for cvx_constant case
%

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

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
