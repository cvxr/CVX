function y = cvx_isconvex( x )
y = false;
for k = 1 : prod( size( x ) ),
    if ~cvx_isconvex( x{k} ), return; end
end
y = true;

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
