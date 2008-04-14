function cvx_optval = sum_square( x, varargin )

%SUM_SQUARE_POS   Internal cvx version.

error( nargchk( 1, 2, nargin ) );
cvx_begin
    variable x2( size( x ) );
    minimize( sum_square( x2, varargin{:} ) );
    x2 >= x;
cvx_end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
