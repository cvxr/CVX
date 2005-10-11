function spy( prob )

global cvx___
prob = index( prob );
p = cvx___.problems( prob );
nn = length( p.reserved );
A = cvx_basis( p.equalities );
if size( A, 2 ) < nn, 
    A( :, nn ) = 0;
end
if ~isempty( p.objective ),
    cc = cvx_basis( p.objective );
    if size( cc, 2 ) < nn, cc( :, nn ) = 0; end
    A = [ -cc ; A ];
end

spy( A );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
