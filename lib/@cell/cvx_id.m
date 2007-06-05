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
