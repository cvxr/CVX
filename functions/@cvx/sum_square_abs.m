function cvx_optval = sum_square_abs( x, varargin )

%SUM_SQUARE_ABS   Internal cvx version.

error( nargchk( 1, 2, nargin ) ); %#ok
cvx_optval = quad_over_lin( x, 1, varargin{:} );

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
