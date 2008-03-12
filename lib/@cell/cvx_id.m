function y = cvx_id( x )
global cvx___
if cvx___.mversion < 7.1,
    nx = prod( size( x ) );
    y = zeros( 1, nx );
    for k = 1 : nx,
        y(k) = cvx_id( x{k} );
    end
else
    y = cellfun( @cvx_id, x );
end
if isempty( y ),
    y = -Inf;
else
    y = max( y( : ) );
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
