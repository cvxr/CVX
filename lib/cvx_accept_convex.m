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

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
