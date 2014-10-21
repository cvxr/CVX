function y = cvx_id( x )
if isempty(x), y = -Inf; return; end
y = max( cellfun( @cvx_id, x ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
