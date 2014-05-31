function v = subsref( x, S )
try
    v = subsref( x.value_, S );
catch
    cvx_throw( 'Invalid tuple reference: %s%s', inputname(1), cvx_subs2str( S ) );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
