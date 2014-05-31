function cvx_optval = norms_largest( varargin )

%NORMS_LARGEST Computation of multiple norm_largest() norms.
%   NORMS_LARGEST( X, K, DIM ) provides a means to compute the largest-k
%   norms of multiple vectors packed into a matrix or N-D vector. This is
%   useful for performing max-of-norms or sum-of-norms calculations.
%
%   If DIM is omitted, the norms are computed along the first non-singleton
%   dimension. 
%
%   See NORM_LARGEST.
%
%   Disciplined convex programming information:
%       NORMS_LARGEST is convex and non-monotonic, so its input must be affine.

[ sx, x, k, dim ] = cvx_get_dimension( varargin, 3 );

if ~isnumeric( k ) || ~isreal( k ) || length( k ) ~= 1,
    cvx_throw( 'Second argument must be a real scalar.' );
end

cvx_optval = sum_largest( abs( x ), k, dim );

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
