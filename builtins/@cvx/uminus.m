function z = uminus( x )

%   Disciplined convex programming information for UMINUS (-):
%      Unary minus may be used in DCPs without restrictions---with the
%      understanding, of course, that it produces a result with the 
%      opposite cuvature to its argument; i.e., the negative of a convex
%      expression is concave, and vice versa.
%
%   Disciplined geometric programming information for UMINUS(-):
%      Negation of non-constant values may not be used in disciplined
%      geometric programs.

z = cvx( x.size_, -x.basis_ );
tt = cvx_isvalid( z, true );
if ~all( tt ),
    cvx_dcp_error( '-', 'unary', cvx_subsref( x, ~tt ) );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
