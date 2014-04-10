function z = mpower( x, y )

%   Disciplined convex/geometric programming information for MPOWER (^):
%      The CVX version of the matrix power function Z=X.^Y supports only
%      the case where X and Y are scalars. In such instances, the rules
%      are identical to those outlined in the help for CVX/POWER.

if length( x ) > 1 && length( y ) > 1,
    error( 'Disciplined convex programming error:\n    Matrix powers not permitted.', 1 ); %#ok
end
try
    z = power( x, y, '^' );
catch exc
	if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
	else rethrow( exc ); end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
