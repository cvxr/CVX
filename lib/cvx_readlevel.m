function y = cvx_readlevel( x )
if ndims( x ) <= 2,
    y = sparse( size( x, 1 ), size( x, 2 ) );
else
    y = zeros( size( x ) );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
