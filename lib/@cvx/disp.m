function disp( x, prefix )
error( cvx_verify( x ) );
if nargin < 2, 
    prefix = ''; 
end
disp( [ prefix, 'cvx affine expression (', type( x ), ')' ] );
dual = getdual( x );
if ~isempty( dual ),
    disp( [ prefix, '   tied to dual variable: ', dual ] );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
