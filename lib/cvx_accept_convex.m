function x = cvx_accept_convex( x )
if isa( x, 'cvx' ),
    t = cvx_vexity( x ) > 0;
    if any( t( : ) ),
        global cvx___
        prob = cvx___.problems(end).self;
        if all( t ),
            src = x;
            dst = newtemp( prob, size( src ) );
            x = dst;
        else
            src = x( t );
            dst = newtemp( prob, size( src ) );
            x( t ) = dst;
        end
        newcnstr( prob, src(:), dst(:), '<=' );
    end
end
