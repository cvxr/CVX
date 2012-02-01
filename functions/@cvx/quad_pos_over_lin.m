function cvx_optval = quad_pos_over_lin( x, y, varargin ) %#ok

%QUAD_POS_OVER_LIN   Internal cvx version.

error( nargchk( 2, 3, nargin ) );
if ~isreal( x ), 
    error( 'First input must be real.' ); 
end
cvx_begin
    variable x2( size(x) )
    minimize quad_over_lin( x2, y, varargin{:} );
    x2 >= x; %#ok
cvx_end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

