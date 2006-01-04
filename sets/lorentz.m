function cvx_optpnt = lorentz( sx, dim, iscplx )
error( nargchk( 1, 3, nargin ) );

%
% Check size vector
%

[ temp, sx ] = cvx_check_dimlist( sx, false );
if ~temp,
    error( 'First argument must be a non-empty dimension vector.' );
end
nd = length( sx );
sy = sx;

%
% Check dimension
%

if nargin < 2 | isempty( dim ),
    dim = [ find( sx > 1 ), 1 ];
    dim = dim( 1 );
elseif ~isnumeric( dim ) | dim < 0 | dim ~= floor( dim ),
    error( 'Second argument must be a dimension (or zero).' );
end
if dim <= 0 | dim > nd | sx( dim ) == 1,
    nv  = 1;
    dim = 0;
else,
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
else,
    iscplx = logical( iscplx );
end

%
% Build the cones
%

nx = prod( sx );
ny = nx / nv;
if iscplx, 
    nv = nv * 2;
    nx = nx * 2;
end
cvx_begin_set
    variables x( nv, ny ) y( 1, ny )
    global cvx___
    p = index( cvx_problem );
    cvx___.problems( p ).reserved( 2 : end ) = 1;
    cvx___.problems( p ).cones = struct( 'type', 'lorentz', 'indices', ...
        [ reshape( 2 : nx + 1, nv, ny ) ; nx + 2 : nx + ny + 1 ] );
cvx_end_set

%
% Permute and reshape as needed
%

if iscplx,
    cvx_optpnt.x = cvx_r2c( cvx_optpnt.x, 1 );
end
nleft = prod( sx( 1 : dim - 1 ) );
if nleft > 1,
    cvx_optpnt.x = reshape( cvx_optpnt.x, [ nv, nleft, ny / nleft ] );
    cvx_optpnt.y = reshape( cvx_optpnt.y, [ 1,  nleft, ny / nleft ] );
    cvx_optpnt.x = permute( cvx_optpnt.x, [ 2, 1, 3 ] );
    cvx_optpnt.y = permute( cvx_optpnt.y, [ 2, 1, 3 ] );
end
cvx_optpnt.x = reshape( cvx_optpnt.x, sx );
cvx_optpnt.y = reshape( cvx_optpnt.y, sy );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
