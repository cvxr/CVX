function cvx_optval = sum_square_abs( x, varargin )

%SUM_SQUARE_ABS   Internal cvx version.

error( nargchk( 1, 2, nargin ) ); %#ok
cvx_optval = quad_over_lin( x, 1, varargin{:} );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
