function x = direction( x )
global cvx___
x = cvx___.problems( index( x ) ).direction;

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
