function minimize( varargin )

%MINIMIZE Specifiies a convex (or affine) objective to be maximized.

if nargin < 1,
    cvx_throw( 'Objective expression missing.' );
elseif iscellstr( varargin ),
    x = evalin( 'caller', sprintf( '%s ', varargin{:} ) );
elseif nargin > 1,
    cvx_throw( 'Too many input arguments.' );
else
    x = varargin{1};
end
evalin( 'caller', 'cvx_verify' );
cvx_pushobj( 'minimize', x );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
