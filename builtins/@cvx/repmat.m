function y = repmat( x, varargin )

%   Disciplined convex/geometric programming information for REPMAT:
%      REPMAT(A,DIMS) imposes no convexity restrictions on A. The
%      dimension arguments must be constant.

s = x.size_;
nx = reshape( 1 : prod( s ), s );
try
    nx = repmat( nx, varargin{:} );
catch
    throw( exc );
end
y = cvx( size( nx ), x.basis_( :, nx ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
