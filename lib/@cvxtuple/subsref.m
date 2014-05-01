function v = subsref( x, S )
try
    v = subsref( x.value_, S );
catch exc    
    tmp = cvx_subs2str( S );
    error( 'CVX:Tuple', 'Invalid tuple reference: %s%s', inputname(1), tmp );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
