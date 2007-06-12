function y = cvx_constant( x )
y = cell( size( x ) );
for k = 1 : prod( size( y ) ),
    y{ k } = cvx_constant( x{ k } );
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
