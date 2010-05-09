function cvx_optval = sum_square( x, varargin )

%SUM_SQUARE   Internal cvx version.

error( nargchk( 1, 2, nargin ) );
if ~isreal( x ),
    error( 'Disciplined convex programming error:\n   The argument to SUM_SQUARE must be real and affine.', 1 ); %#ok
end
cvx_optval = quad_over_lin( x, 1, varargin{:} );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
