function z = iscurrent( prob )
global cvx___
if isempty( cvx___.stack ),
    z = false;
else,
    z = cvx___.stack{ end } == prob;
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
