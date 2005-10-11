function cvx_optval = abs( x )
    
sx = size( x );
cvx_begin
    variables y( sx );
    minimize y;
    { x, y } == lorentz( sx, 0, ~isreal( x ) );
cvx_end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
