function cvx_optpnt = lorentz( sx, dim, iscplx )
error( nargchk( 1, 3, nargin ) );

%
% Check size vector
%

[ temp, sx ] = cvx_check_dimlist( sx, true );
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
end
sy = sx;
nd = length( sx );
if dim <= 0 | dim > nd | sx( dim ) == 1,
    nv  = 1;
    dim = 0;
else
    nv = sx( dim );
    sy( dim ) = 1;
end

%
% Check complex flag
%

if nargin < 3 | isempty( iscplx ),
    iscplx = false;
elseif length( iscplx ) ~= 1,
    error( 'Third argument must be a scalar.' );
else
    iscplx = logical( iscplx );
end

%
% Quick exit for the length-0 case
%

if nv == 0,
    cvx_optpnt.y = cvx_collapse( nonnegative( sy ), false );
    cvx_optpnt.x = zeros( sx );
    return
end

%
% Build the cone
%

nx = prod( sx );
ny = nx / nv;
if iscplx,
    nv = nv * 2;
    nx = nx * 2;
end
cvx_begin_set
    variables x( nv, ny ) y( 1, ny )
    [ tx, dummy ] = find( cvx_basis( x ) );
    [ ty, dummy ] = find( cvx_basis( y ) );
    newnonl( cvx_problem, 'lorentz', [ reshape( tx, nv, ny ) ; reshape( ty, 1, ny ) ] );
cvx_end_set

%
% Permute and reshape as needed
%

if iscplx,
    x = cvx_r2c( x, 1 );
end
nleft = prod( sx( 1 : dim - 1 ) );
if nleft > 1,
    x = reshape( x, [ nv, nleft, ny / nleft ] );
    y = reshape( y, [ 1,  nleft, ny / nleft ] );
    x = permute( x, [ 2, 1, 3 ] );
    y = permute( y, [ 2, 1, 3 ] );
end
x = reshape( x, sx );
y = reshape( y, sy );
cvx_optpnt = cvxtuple( struct( 'x', x, 'y', y ) );

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
