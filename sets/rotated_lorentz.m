function set = rotated_lorentz( varargin )

%ROTATED_LORENTZ   Rotated real second-order cone.
%   ROTATED_LORENTZ(N), where N is a positive integer, creates a column
%   variable of length N and two scalar variables, and constrains them
%   to lie in a rotated second-order cone. That is, given the declaration
%       variables x(n) y z
%   the constraint
%       {x,y,z} == rotated_lorentz(n)
%   is equivalent to
%       norm(x,2) <= 2*geo_mean([y,z])
%   except that using ROTATED_LORENTZ is more efficient.
%
%   ROTATED_LORENTZ(SX,DIM), where SX is a valid size vector and DIM is a
%   positive integer, creates an array variable of size SX and two 
%   array variables of size SY (see below) and applies the second-order cone
%   constraint along dimension DIM. That is, given the declarations
%       sy = sx; sy(min(dim,length(sx)+1))=1;
%       variables x(sx) y(sy) z(sz)
%   the constraint
%       {x,y,z} == rotated_lorentz(sx,dim)
%   is equivalent to
%       norms(x,2,dim) <= 2*geo_mean(cat(dim,y,z),dim)
%   except, again, ROTATED_LORENTZ is more efficient. DIM is optional; if
%   it is omitted, the first non-singleton dimension is used.
%
%   ROTATED_LORENTZ(SX,DIM,CPLX) creates real second-order cones if CPLX is
%   FALSE, and complex second-order cones if CPLX is TRUE. The latter case
%   is equivalent to COMPLEX_ROTATED_LORENTZ(SX,DIM).
%
%   Disciplined convex programming information:
%       ROTATED_LORENTZ is a cvx set specification. See the user guide for
%       details on how to use sets.

[ sx, dim, iscplx ] = cvx_get_dimension( varargin, 2, 'nox', true, 'zero', true );
if isempty( iscplx ),
    iscplx = false;
elseif length( iscplx ) ~= 1,
    cvx_throw( 'Third argument must be a scalar.' );
else
    iscplx = logical( iscplx );
end
sy = sx;
nv = sx( dim );
sy( dim ) = 1;
ny = prod( sy );
if nv == 0 || ny == 0,
    x = cvx( sx, [] );
    if ny == 0,
        y = cvx( sy, [] );
        z = cvx( sy, [] );
    else
        cone = nonnegative( 2 * prod(sy) );
        y = reshape( cone(1:2:end), sy );
        z = reshape( cone(2:2:end), sy );
    end
else
    if iscplx, 
        nv = nv * 2; 
    end
    cvx_begin set
        variable x( nv, ny ) 
        variable y( 1,  ny ) nonnegative_
        variable z( 1,  ny ) nonnegative_
        cvx_pushcone( true, 'rotated_lorentz', [ x ; y ; z ] ); %#ok
    cvx_end
    if iscplx,
        x = x(1:2:end,:) + 1j * x(2:2:end,:);
        nv = nv * 0.5;
    end
    if sx( dim ) > 1,
        nleft = prod( sx( 1 : dim - 1 ) );
        if nleft > 1,
            x = reshape( x, [ nv, nleft, ny / nleft ] );
            x = permute( x, [ 2, 1, 3 ] );
        end
    end
    x = reshape( x, sx );
    y = reshape( y, sy );
    z = reshape( z, sy );
end
set = cvxtuple( struct( 'x', x, 'y', y, 'z', z ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
