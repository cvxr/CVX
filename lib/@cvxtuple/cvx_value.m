function y = cvx_value( x )
y = do_value( x.value_ );
function y = do_value( x )
switch class( x ),
    case 'struct', y = cell2struct( do_value( struct2cell( x ) ), fieldnames( x ), 1 );
    case 'cell',   y = cellfun( @do_value, x, 'UniformOutput', false );
    otherwise,     y = cvx_value( x );
end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
