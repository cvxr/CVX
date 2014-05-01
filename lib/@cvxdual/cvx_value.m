function v = cvx_value( x )
global cvx___
try
    tv = subsref( cvx___.problems( x.problem_ ).duals, x.name_ );
catch
    v = NaN;
    return
end
v = value_( tv, cvx___.y );
function y = value_( x, y )
switch class( x ),
    case 'struct', y = cell2struct( value_( struct2cell( x ), y ), fieldnames( x ), 1 );
    case 'cell',   y = cellfun( @(z)value_( z, y ), x, 'UniformOutput', false );
    case 'cvx',    y = value( x, y );
    otherwise,     y = x;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
