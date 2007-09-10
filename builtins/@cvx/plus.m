function z = plus( x, y, isdiff, cheat )

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

%
% Default arguments
%

if nargin < 4,
    if nargin < 3,
        isdiff = false;
    end
    cheat = false;
end

%
% Check sizes
%

sx = size( x );
sy = size( y );
xs = all( sx == 1 );
ys = all( sy == 1 );
if xs,
    sz = sy;
elseif ys,
    sz = sx;
elseif ~isequal( sx, sy ),
    error( 'Matrix dimensions must agree.' );
else
    sz = sx;
end

%
% Check vexity
%

if ~cheat,
    temp = cvx_vexity( x ) .* cvx_vexity( y );
    if isdiff,
        bad = temp( : ) > 0;
    else
        bad = temp( : ) < 0;
    end
    if any( bad ),
        if isdiff,
            error( sprintf( 'Disciplined convex programming error:\n   Differences between convex or concave terms are forbidden.' ) );
        else
            error( sprintf( 'Disciplined convex programming error:\n   Sums of convex and concave terms are forbidden.' ) );
        end
    end
end

%
% Apply operation, stretching basis matrices as needed
%

if any( sz == 0 ),
    bz = sparse( 1, 0 );
else
    x  = cvx( x );
    y  = cvx( y );
    bx = x.basis_;
    by = y.basis_;
    if isdiff,
        by = -by;
    end
    if xs & ~ys,
        nz = prod( sz );
        bx = bx( :, ones( 1, nz ) );
    elseif ys & ~xs,
        nz = prod( sz );
        by = by( :, ones( 1, nz ) );
    end
    nx = size( bx, 1 );
    ny = size( by, 1 );
    if nx < ny,
        if issparse( by ) & ~issparse( bx ), bx = sparse( bx ); end
        bx( ny, : ) = 0;
    elseif ny < nx,
        if issparse( bx ) & ~issparse( by ), by = sparse( by ); end
        by( nx, : ) = 0;
    end
    bz = bx + by;
end

%
% Construct result
%

z = cvx( sz, bz );

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
