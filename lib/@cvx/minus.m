function z = minus( x, y, cheat )
if nargin < 3, cheat = false; end
z = plus( x, y, true, cheat );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

