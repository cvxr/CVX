function s = cvx_verify( x )

if valid( x ),
    s = '';
else,
    s = 'Attempt to reference a completed cvx specification.';
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

