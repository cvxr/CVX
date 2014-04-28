function disp( x, prefix )
if nargin < 2,
    prefix = '';
end
str = [ prefix, 'cvx: ', cvx_class( x, true, true, true ) ];
if numel( x ) > 1,
    str = [ str, ' ', type( x ) ];
elseif any( strcmp( str(end-5:end), { 'nstant', ': zero' } ) ),
    str = [ prefix, 'cvx: ', num2str(cvx_constant(x)) ];
end
disp( str );
dual = cvx_getdual( x );
if ~isempty( dual ) && isstruct( dual ),
    disp( [ prefix, '   tied to dual variable: ', dual.subs ] );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
