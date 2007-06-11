function cvx_optval = quad_pos_over_lin( x, y, dim )

%QUAD_POS_OVER_LIN   Internal cvx version.

error( nargchk( 2, 3, nargin ) );
if ~isreal( x ), 
    error( 'First input must be real.' ); 
end

sx = size( x );
if nargin < 3,
    dim = cvx_default_dimension( sx );
end

cvx_begin
    variable x2( sx )
    minimize quad_over_lin( x2, y, dim );
    x2 >= x;
cvx_end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

