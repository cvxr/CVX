function y = cvx_default_dimension( sx )

% CVX_DEFAULT_DIMENSION   Default dimension for SUM, MAX, etc. 

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

y = find( sx ~= 1 );
if isempty( y ), 
    y = 1; 
else
    y = y( 1 ); 
end


