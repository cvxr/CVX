function y = cvx_id( x )
global cvx___
if cvx___.hcellfun,
    y = cellfun( @cvx_id, x );
else
    nx = prod( size( x ) );
    y = zeros( 1, nx );
    for k = 1 : nx,
        y(k) = cvx_id( x{k} );
    end
end
if isempty( y ),
    y = -Inf;
else
    y = max( y( : ) );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
