function z = norm_nuc( x )

%NORM_NUC   Nuclear norm of a matrix.
%   NORM_NUC(X) = SUM(SVD(X)). X must be a 2-D matrix, real or complex.
%
%   Disciplined convex programming information:
%       NORM_NUC(X) is convex and nonmontonic in X, so X must be affine.

error( nargchk( 1, 1, nargin ) );
z = sum(svd(x));

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
