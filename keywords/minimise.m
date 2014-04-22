function minimise( varargin )

%MINIMISE Specifiies a convex (or affine) objective to be maximized.

prob = evalin( 'caller', 'cvx_verify' );
if nargin < 1,
    error( 'Objective expression missing.' );
elseif iscellstr( varargin ),
    x = evalin( 'caller', sprintf( '%s ', varargin{:} ) );
elseif nargin > 1,
    error( 'Too many input arguments.' );
else
    x = varargin{1};
end
try
    newobj( prob, 'minimize', x );
catch exc
    rethrow( exc )
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
