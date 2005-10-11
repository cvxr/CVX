function y = permute( x, order )
error( cvx_verify( x ) );
sx = size( x );
try,
    ndxs = builtin( 'permute', reshape( 1 : prod( sx ), sx ), order );
catch,
    error( lasterr );
end
y = reshape( cvx_subref( x, ndxs( : ) ), size( ndxs ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
