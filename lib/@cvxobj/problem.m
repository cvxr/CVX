function y = problem( x )
global cvx___
y = cvx___.problems( x.index_ ).self;
if x.index_ > length( cvx___.problems ) | x.id_ < id( y ),
    error( 'Object references a cvx specification that has already been completed.' );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
