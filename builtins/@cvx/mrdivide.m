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

z = mtimes( x, y, 'rdivide' );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
