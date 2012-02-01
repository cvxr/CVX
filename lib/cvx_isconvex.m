function y = cvx_isconvex( x, full )
error( nargchk( 1, 2, nargin ) );
if nargin == 2,
    y = ~imag( x );
elseif isreal( x ),
    y = true;
else
    y = nnz(imag(x)) == 0;
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
