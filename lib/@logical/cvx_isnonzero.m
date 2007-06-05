function y = cvx_isnonzero( x, full )
error( nargchk( 1, 2, nargin ) );
y = x( : ) ~= 0;
if nargin < 2,
    y = all( y );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
