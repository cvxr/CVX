function ans = cvx_isvalid( x )
ans = false;
for k = 1 : prod( size( x ) ),
    if ~cvx_isvalid( x{k} ), 
        return; 
    end
end
ans = true;

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
