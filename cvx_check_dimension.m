function y = cvx_check_dimension( x, zero_ok )

% CVX_CHECK_DIMENSION	Verifies that the input is valid dimension.

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

if isnumeric( x ) & length( x ) == 1 & isreal( x ) & x < Inf & x == floor( x ),
    y = x > 0 | zero_ok;
else,
    y = 0;
end
