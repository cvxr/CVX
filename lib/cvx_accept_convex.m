function x = cvx_accept_convex( x )
t = cvx_vexity( x ) > 0;
if any( t( : ) ),
    prob = cvxprob( 'current' );
    if all( t ),
        src = x;
        dst = newtemp( prob, size( src ) );
        x = dst;
    else,
        src = x( t );
        dst = newtemp( prob, size( src ) );
        x( t ) = dst;
    end
    newcnstr( prob, src(:), dst(:), '<=' );
end
