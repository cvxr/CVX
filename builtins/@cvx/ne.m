function a = ne( x, y )

%Disciplined convex/geometric programming information for NE (~=):
%   Not-equal constraints violate both the DCP and DGP rulesets. Thus 
%   not-equal expressions may only appear in CVX models when both sides
%   are constant.

error( nargchk( 2, 2, nargin ) );
try
    newcnstr( evalin( 'caller', 'cvx_problem', '[]' ), x, y, '~=' );
catch
    error( cvx_lasterr );
end
if nargout > 0,
    a = 'Constraint accepted';
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
