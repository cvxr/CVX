function z = plus( x, y )

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

persistent params
if isempty( params ),
    params.map   = {};
    params.funcs = { @plus_nc };
    params.name  = '+';
end

try
    z = cvx_binary_op( params, x, y );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
