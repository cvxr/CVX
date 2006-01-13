function cvx_optpnt = rotated_lorentz( sx, dim, iscplx )
error( nargchk( 1, 3, nargin ) );

%
% Check size vector
%

[ temp, sx ] = cvx_check_dimlist( sx, false );
if ~temp,
    error( 'First argument must be a non-empty dimension vector.' );
end

%
% Check dimension
%

if nargin < 2 | isempty( dim ),
    dim = cvx_default_dimension( sx );
elseif ~cvx_check_dimension( dim, true ),
    error( 'Second argument must be a dimension (or zero).' );
elseif dim == 0 | dim > nd | sx( dim ) == 1,
    dim = cvx_default_dimension( sx );
end
nd = length( sx );
if dim > nd,
    sx( end + 1 : dim ) = 1;
    nd = dim;
end

%
% Check complex flag
%

if nargin < 3 | isempty( iscplx ),
    iscplx = false;
elseif length( iscplx ) ~= 1,
    error( 'Third argument must be a scalar.' );
else,
    iscplx = logical( iscplx );
end

%
% Build the cvx module
%

if iscplx,
    sx( dim ) = 2 * sx( dim );
end

if sx( dim ) == 1,
    cone = semidefinite( [ 2, 2, sx ] );
    cvx_optpnt.x = reshape( cone( 2, 1, : ), sx );
    cvx_optpnt.y = reshape( cone( 1, 1, : ), sx );
    cvx_optpnt.z = reshape( cone( 2, 2, : ), sx );
else,
    sx( dim ) = sx( dim ) + 1;
    cone = lorentz( sx, dim );
    ndxs = cell( 1, nd );
    [ ndxs{:} ] = deal( ':' );
    ndxs{ dim } = 1 : sx( dim ) - 1;
    cvx_optpnt.x = cone.x( ndxs{:} );
    ndxs{ dim } = sx( dim );
    temp = cone.x( ndxs{:} );
    cvx_optpnt.y = cone.y + temp;
    cvx_optpnt.z = cone.y - temp;
end

if iscplx,
    cvx_optpnt.x = cvx_r2c( cvx_optpnt.x, dim );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
