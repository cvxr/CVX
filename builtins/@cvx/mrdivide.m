function z = mrdivide( x, y )

%   Disciplined convex/geomtric programming information for MRDIVIDE:
%      The MRDIVIDE operation X/Y is quite often employed with Y as a 
%      scalar. In that case, it is equivalent to the RDIVIDE operation
%      X./Y, and must obey the same rules as outlined in the help for 
%      CVX/RDIVIDE.
%   
%      When Y is a matrix, the MRDIVIDE operation X/Y is equivalent to
%      X*inv(Y) for both DCP and DGP purposes. The inv() operation is 
%      not supported for non-constant expressions, so Y must be both 
%      constant and nonsingular. The resulting matrix multiplication 
%      must obey the same rules as outlined in the help for CVX/MTIMES.

try
    sz = size( y );
    if all( sz == 1 ),
        z = rdivide( x, y, '/' );
    elseif length( sz ) > 2,
        error( 'Inputs must be 2-D, or at least one input must be scalar.' );
    elseif sz( 1 ) ~= sz( 2 ) && length( sz ) == 2,
        error( 'Non-square matrix divisors are not supported in CVX.' );
    elseif ~cvx_isconstant( y ),
        error( 'Matrix divisors must be constant.' );
    else
        z = mtimes( x, cvx_constant( y ), '/' );
    end
catch exc
	if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
	else rethrow( exc ); end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
