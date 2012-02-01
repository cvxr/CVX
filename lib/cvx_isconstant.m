function y = cvx_isconstant( x, full )
error( nargchk( 1, 2, nargin ) );
if nargin == 1,
    y = true;
else
    y = true( size( x ) );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
