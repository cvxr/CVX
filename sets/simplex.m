function cvx_optpnt = simplex( sx, dim )
error( nargchk( 1, 2, nargin ) );

%
% Check size vector
%

[ temp, sx ] = cvx_check_dimlist( sx, false );
if ~temp,
    error( 'First argument must be a non-empty dimension vector.' );
end
nd = length( sx );

%
% Check dimension
%

if nargin < 2 | isempty( dim ),
    dim = [ find( sx > 1 ), 1 ];
    dim = dim( 1 );
elseif ~isnumeric( dim ) | dim < 0 | dim ~= floor( dim ),
    error( 'Second argument must be a dimension.' );
elseif dim > nd,
    sx( end + 1 : dim ) = 1;
    nd = dim;
end

%
% Construct set
%

cvx_begin_set
   variables x( sx )
   sum( x, dim ) == 1
   x >= 0
cvx_end_set

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
