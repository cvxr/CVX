function cvx_optval = square_pos( x )

%SQUARE_POS   Internal cvx version.

error( nargchk( 1, 1, nargin ) );
cvx_optval = pow_pos( x, 2 );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
