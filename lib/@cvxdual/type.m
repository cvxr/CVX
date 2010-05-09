function st = type( x )

s   = size( x );
len = prod( s );
isr = isreal( x );

if len == 0,
    st = 'unassigned';
elseif len == 1,
    isstruct = 0;
    if isr,
        st = 'scalar';
    else
        st = 'complex scalar';
    end
else
    dof = size( cvx_basis( x ), 2 ) - 1;
    isstruct = dof < ( 2 - isr ) * len;
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
    if isstruct,
        st = sprintf( '%s (%d d.o.f.)', st, dof );
    end
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
