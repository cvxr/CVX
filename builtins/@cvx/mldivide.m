function z = mldivide( x, y )

%   Disciplined convex/geomtric programming information for MLDIVIDE:
%      The MLDIVIDE operation X\Y can be employed with X as a scalar. In
%      that case, it is equivalent to the LDIVIDE operation X.\Y, and 
%      must obey the same rules as outlined in the help for CVX/LDIVIDE.
%   
%      When X is a matrix, the MRDIVIDE operation X\Y is equivalent to
%      inv(X)*Y for both DCP and DGP purposes. The inv() operation is 
%      not supported for non-constant expressions, so X must be both 
%      constant and nonsingular. The resulting matrix multiplication 
%      must obey the same rules as outlined in the help for CVX/MTIMES.

sz = size( x );
if all( sz == 1 ),
    z = ldivide( x, y, '\' );
elseif length( sz ) > 2,
    cvx_throw( 'Inputs must be 2-D, or at least one input must be scalar.' );
elseif sz( 1 ) ~= sz( 2 ),
    cvx_throw( 'Non-square matrix divisors are not supported in CVX.' );
elseif ~cvx_isconstant( x ),
    cvx_throw( 'Matrix divisors must be constant.' );
else
    z = mtimes( cvx_constant( x ), y, '\' );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
