function z = plus( x, y, op, cheat )

%   Disciplined convex programming information for PLUS:
%      Both terms in a sum must have the same curvature. Real affine
%      expressions are both convex and concave, so they can be added to
%      any nonlinear expressions. Complex affine (or constant)
%      expressions, however, can only be added to other affine 
%      expressions. So, for example, the following sums are valid:
%         {convex}+{convex}   {concave}+{concave}   {affine}+{affine}
%      The following are not:
%         {convex}+{concave}  {convex}+{complex constant}
%   
%   Disciplined geometric programming information for PLUS:
%      Only log-convex terms may be summed; this includes positive 
%      constants, monomials, posynomials, and generalized posynomials.
%   
%   For vectors, matrices, and arrays, these rules are verified 
%   indepdently for each element.

% Determine sizes

try
    [ x, y, sz, xs, ys ] = cvx_broadcast( x, y );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end
nn = prod( sz );
if ~nn, 
    z = zeros( sz ); 
    return; 
end

% Build basis

x  = cvx( x );
y  = cvx( y );
bx = x.basis_;
by = y.basis_;
if xs,
    bx = bx( :, ones( 1, nn ) );
elseif ys,
    by = by( :, ones( 1, nn ) );
end
[ nx, nv ] = size( bx );
ny = size( by, 1 );
if nx < ny
    if issparse( by ), bx = sparse( bx ); end
    bx = [ bx ; sparse( ny - nx, nv ) ];
elseif ny < nx,
    if issparse( bx ), by = sparse( by ); end
    by = [ by ; sparse( nx - ny, nv ) ];
end
if nargin >= 3 && op == '-',
    bz = bx - by;
else
    bz = bx + by;
end

% Build object

z = cvx( sz, bz );

% Check

if nargin < 4 || ~cheat,
    if nargin < 3, 
        op = '+'; 
    end
    try
        cvx_dcp_error( x, y, z, op );
    catch exc
        throw( exc );
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
