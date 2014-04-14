function st = type( x )

x = cvxaff( x );
if isa( x, 'double' ),
    st = 'unassigned';
    return
end
if iscell( x ),
    strs = cell(1,numel(x));
    for k = 1 : numel(x),
        strs{k} = type(x{k});
    end
    strs = sprintf( '%s, ', strs{:} );
    st = strs(1:end-2);
    return
end
s   = size( x );
len = prod( s );
isr = isreal( x );
if len == 1,
    if isr,
        st = 'scalar';
    else
        st = 'complex scalar';
    end
else
    st = sprintf( '%dx', s );
    st = st( 1 : end - 1 );
    if ~isr,
        st = [ st, ' complex' ];
    end
    if length( s ) > 2,
        st = [ st, ' array' ];
    elseif any( s == 1 ),
        st = [ st, ' vector' ];
    else
        st = [ st, ' matrix' ];
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
