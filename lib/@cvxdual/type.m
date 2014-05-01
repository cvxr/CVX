function st = type( x )
global cvx___
try
    v = subsref( cvx___.problems( x.problem_ ).duals, x.name_ );
catch
    st = 'invalid';
    return
end
switch class(v),
    case 'cvx',
        st = type( v );
    case { 'cell', 'struct' },
        st = 'composite';
    otherwise,
        st = 'unassigned';
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
