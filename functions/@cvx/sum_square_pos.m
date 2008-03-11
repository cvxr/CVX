function cvx_optval = sum_square( x, dim )

%SUM_SQUARE_POS   Internal cvx version.

error( nargchk( 1, 2, nargin ) );
sx = size( x );
if nargin < 2,
    dim = cvx_default_dimension( sx );
end

cvx_begin
    variable x2( sx )
    minimize sum_square( x2, dim )
    x2 >= x;
cvx_end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
