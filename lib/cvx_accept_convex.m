function x = cvx_accept_convex( x )
global cvx___
if isa( x, 'cvx' ),
    t = cvx_vexity( x ) > 0;
    if any( t( : ) ),
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

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
