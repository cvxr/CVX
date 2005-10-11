function st = type( x, mod )
error( cvx_verify( x ) );

df  = dof( x );
s   = size( x );
len = prod( s );
isr = isreal( x );

if len == 1,
    isstruct = 0;
    if ~isr,
        st = 'complex scalar';
    else,
        st = 'scalar';
    end
else,
    isstruct = df < ( 2 - isr ) * len;
    nd = length( s );
    st = sprintf( '%dx', s );
    st = st( 1 : end - 1 );
    if ~isreal( x ),
        st = [ st, ' complex' ];
    end
    if nd > 2,
    	st = [ st, ' array' ];
    elseif any( s == 1 ),
    	st = [ st, ' vector' ];
    else,
    	st = [ st, ' matrix' ];
    end
end

if cvx_isconstant( x ),
    st = [ st, ' constant' ];
end

if isstruct,
    st = sprintf( '%s: %d d.o.f.', st, df );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
