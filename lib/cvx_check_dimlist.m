function [ y, x ] = cvx_check_dimlist( x, emptyok )

% CVX_CHECK_DIMLIST Verifies the input is a valid dimension list.

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

if nargin < 2 || emptyok,
    xmin = 0;
else
    xmin = 1;
end
if isa( x, 'cell' ),
    nel = numel( x );
    xnew = zeros( 1, nel );
    y = false;
    for k = 1 : nel,
        if ~isnumeric( x{k} ) || length( x{k} ) ~= 1, return; end
        xnew( k ) = x{k};
    end
    x = xnew;
end
if isnumeric( x ) && ndims( x ) <= 2 && ~all( size( x ) > 1 ) && isreal( x ) && ~any( isnan( x ) | isinf( x ) ) && all( x >= xmin ) && all( floor( x ) == x ),
    y = true;
else
    y = false;
end
if y && nargout > 1,
    x = [ x( : )', 1, 1 ];
    x = x( 1 : max( [ 2, find( x > 1 ) ] ) );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
