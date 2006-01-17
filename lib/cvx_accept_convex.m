function x = cvx_accept_convex( x )
t = cvx_vexity( x ) == +1;
if any( t( : ) ),
    prob = cvxprob( 'current' );
    if all( t ),
        src = x;
        dst = newtemp( prob, x );
        x = dst;
    else,
        src = x( t );
        dst = newtemp( prob, src );
        x( t ) = dst;
    end
    newcnstr( prob, src, dst, '>=' );
end
