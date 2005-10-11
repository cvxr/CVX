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

%
% Check dimension
%

if nargin < 2 | isempty( dim ),
    dim = [ find( sx > 1 ), 1 ];
    dim = dim( 1 );
elseif ~isnumeric( dim ) | dim < 0 | dim ~= floor( dim ),
    error( 'Second argument must be a dimension (or zero).' );
elseif dim > nd,
    sx( end + 1 : dim ) = 1;
    nd = dim;
elseif dim == 0,
    dim = min( find( sx == 1 ) );
    if isempty( dim ), dim = nd + 1; end
    sx( dim ) = 1;
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

sd = sx( dim );
nx = prod( sx );
sv = nx / sd;
sy = sx;
sy( dim ) = 1;
if iscplx, 
    sd = sd * 2;
    nx = nx * 2;
end
use_lp = sd == 1;

if use_lp,
    
    cvx_begin_set
        variables x( sx ) y( sy )
        +x <= y;
        -x <= y;
    cvx_end_set
    
else,

    cvx_begin_set
        if iscplx,
            variable x( sx ) complex
            sx( dim ) = sd;
        else,
            variable x( sx )
        end
        variable y( sy )
        tx = reshape( 2 : nx + 1, sx );
        ty = reshape( nx + 2 : nx + sv + 1, sy );
        if sd > 1 & any( sx( 1 : dim - 1 ) > 1 ),
            perm = [ dim, 1 : dim - 1, dim + 1 : nd ];
            tx = permute( tx, perm );
            ty = permute( ty, perm );
        end
        tx = reshape( tx, sd, sv );
        ty = reshape( ty, 1,  sv );
        global cvx___
        p = index( cvx_problem );
        cvx___.problems( p ).reserved( 2 : end ) = 1;
        cvx___.problems( p ).cones = struct( 'type', 'lorentz', 'indices', [ tx ; ty ] );
    cvx_end_set
    
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
