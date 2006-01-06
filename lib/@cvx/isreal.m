function y = isreal( x, full )
error( cvx_verify( x ) );
y = cvx_basis( x );
if nargin > 1 & full,
    y = any( imag( y ), 2 );
else,
    y = isreal( y );
end
    
% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
