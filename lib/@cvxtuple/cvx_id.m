function y = cvx_id( x )
y = do_id( x.value_ );
function y = do_id( x )
switch class( x ),
    case 'struct', y = max( cellfun( @do_id, struct2cell( x ) ) );
    case 'cell',   y = max( cellfun( @do_id, x ) );
    otherwise,     y = cvx_id( x );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
