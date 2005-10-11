function y = issimple( x )
error( cvx_verify( x ) );
b = cvx_basis( x );
if isreal( x ),
    y = all( sum( b ~= 0, 2 ) <= ( b( :, 1 ) ~= 0 ) + 1 );
else,
    bR = real( b );
    bI = imag( b );
    y = all( sum( bR ~= 0, 2 ) <= ( bR( :, 1 ) ~= 0 ) + 1 ) & all( sum( bI ~= 0, 2 ) <= ( bI( :, 1 ) ~= 0 ) + 1 );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
