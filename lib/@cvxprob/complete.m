function y = iscomplete( x )
global cvx___
y = cvx___.problems( index( x ) ).complete;

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
