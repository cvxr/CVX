function s = cvx_verify( x )
global cvx___
for k = 1 : prod(size(x)),
    s = cvx_verify( x{k} );
    if ~isempty( s ), return; end
end
s = '';

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

