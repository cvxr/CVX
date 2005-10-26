function x = subsasgn( x, S, y )
error( nargchk( 3, 3, nargin ) );
error( cvx_verify( x, y ) );

%
% Test subscripts
%

szx = size( x );
szy = size( y );
nlx = prod( szx );
try,
    temp = reshape( 1 : nlx, szx );
    ndx_x = builtin( 'subsasgn', temp, S, zeros( szy ) );
catch,
    error( lasterr );
end
szx_n = size( ndx_x );

%
% Assign data
%

[ prob, x, y ] = cvx_operate( [], x, y );
bx = x.basis_;
if any( szx_n < szx ),
    bx = bx( ndx_x, : );
else,
    if any( szx_n > szx ),
        bx( end + 1, : ) = 0;
        ndx_x( ndx_x == 0 ) = size( bx, 1 ) + 1;
        bx = bx( ndx_x, : );
        temp = reshape( 1 : prod( szx_n ), szx_n );
    end
    ndx_x = builtin( 'subsref', temp, S );
    ndx_x = ndx_x( : );
    nlz = length( ndx_x );
    by = y.basis_;
    if size( by, 2 ) < nlz,
        by = by( ones( 1, nlz ), : );
    end
    bx( ndx_x, 1 : size( by, 2 ) ) = by;
end

%
% Create the new object
%

x = cvx( prob, szx_n, bx );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
