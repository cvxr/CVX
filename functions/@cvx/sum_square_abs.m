function cvx_optval = sum_square_abs( x, varargin )

%SUM_SQUARE_ABS   Internal cvx version.

error( nargchk( 1, 2, nargin ) );
cvx_optval = quad_over_lin( x, 1, varargin{:} );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
