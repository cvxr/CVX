function x = cvx_accept_convex( x )
t = cvx_isconvex( x, true );
if any( t ),
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
