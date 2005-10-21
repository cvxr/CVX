function st = type( x, mod )
error( cvx_verify( x ) );

df  = x.dof_;
s   = x.size_;
len = prod( s );
isr = isreal( x.basis_ );

if len == 1,
    isstruct = 0;
    if ~isr,
        st = 'complex scalar';
    else,
        st = 'scalar';
    end
else,
    if isempty( df ),
        isstruct = false;
    else,
        isstruct = df < ( 2 - isr ) * len;
    end
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
