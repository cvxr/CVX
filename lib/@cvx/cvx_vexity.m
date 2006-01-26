function v = cvx_vexity( x )
error( cvx_verify( x ) );

global cvx___
p = cvx___.problems( index( problem( x ) ) );
b = cvx_basis( x );
s = size( x );

temp1 = p.vexity( 1 : size( b, 2 ) );
temp2 = full( b * temp1( : ) );
temp3 = b( :, temp1 ~= 0 );
temp3 = full( sum( reshape( abs( temp3 ), size( temp3 ) ), 2 ) );
temp2( abs( temp2 ) ~= temp3 ) = NaN;
v = reshape( sign( temp2 ), s );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
