function z = sigma_max( x )

%SIGMA_MAX    Maximum singular value.
%   SIGMA_MAX(X) returns the maximum singular value of X. X must be a 2-D
%   matrix, real or complex. SIGMA_MAX(X) is synonymous with NORM(X).
%
%   Disciplined convex programming information:
%       SIGMA_MAX(X) is convex and nonmontonic in X, so X must be affine.

error( nargchk( 1, 1, nargin ) );
z = norm( x, 2 );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
