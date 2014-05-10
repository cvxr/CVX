function a = ne( x, y )

%   Disciplined convex/geometric programming information for NE (~=):
%      Not-equal constraints violate both the DCP and DGP rulesets. Thus
%      not-equal expressions may only appear in CVX models when both 
%      sides are constant.

try
	evalin( 'caller', 'cvx_verify' );
	b = cvx_pushcnstr( x, y, '~=' );
	if nargout, a = b; end
catch exc
	if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
	else rethrow( exc ); end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
