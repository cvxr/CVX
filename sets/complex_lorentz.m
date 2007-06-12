function cvx_optpnt = complex_lorentz( sx, dim )
error( nargchk( 1, 2, nargin ) );
if nargin == 1,
    cvx_optpnt = lorentz( sx, [], true );
else
    cvx_optpnt = lorentz( sx, dim, true );
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
