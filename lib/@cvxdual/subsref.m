function v = subsref( x, S )
global cvx___
if numel( S ) ~= 1 || ~isequal( S.type, '{}' ),
    cvx_throw( 'Only cell subscripting allowed for dual variables.' );
end
S = [ x.name_, S ];
try
    subsref( cvx___.problems( x.problem_ ).duals, S );
catch exc
    tmp = cvx_subs2str( S );
    cvx_throw( 'Invalid dual variable: %s', tmp(2:end) );
end
v = cvxdual( x.problem_, S );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
