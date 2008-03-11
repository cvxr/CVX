function st = type( x, usegeo )

df  = x.dof_;
if nargin > 1 & usegeo,
    geo = any( df < 0 );
else
    geo = false;
end
df  = abs( df );
s   = x.size_;
len = prod( s );
isr = isreal( x.basis_ );

if len == 1,
    isstruct = 0;
    nzs = len;
    if ~isr,
        st = 'complex scalar';
    else
        st = 'scalar';
    end
else
    isstruct = ~isempty( df ) & ( sum( df ) < ( 2 - isr ) * len );
    if ~isstruct,
        nzs = nnz( any( x.basis_, 1 ) );
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
    else
        st = [ st, ' matrix' ];
    end
end

%if cvx_isconstant( x ),
%    st = [ st, ' constant' ];
%end
if isstruct,
    st = sprintf( '%s, %d d.o.f.', st, df );
elseif nzs < len,
    if nzs > 1,
        st = sprintf( '%s, %d nonzeros', st, nzs );
    elseif nzs == 1,
        st = sprintf( '%s, 1 nonzero', st );
    end
end
if geo,
    st = [ st, ', geometric' ];
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
