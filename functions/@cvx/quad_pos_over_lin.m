function cvx_optval = quad_pos_over_lin( x, y, varargin )

%QUAD_POS_OVER_LIN   Internal cvx version.

error( nargchk( 2, 3, nargin ) ); %#ok
if ~isreal( x ),
    error( 'Disciplined convex programming error:\n   The argument to QUAD_POS_OVER_LIN must be real.', 1 ); %#ok
end
if all( cvx_sign( x ) ) > 0,
    cvx_optval = quad_over_lin( x, y, varargin{:} );
else
    x2 = [];
    cvx_begin
        variable x2( size(x) )
        minimize quad_over_lin( x2, y, varargin{:} );
        x2 >= x; %#ok
    cvx_end
end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

