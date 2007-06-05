function x = subsasgn( x, S, y )
error( nargchk( 3, 3, nargin ) );

%
% Test subscripts
%

szx = size( x );
szy = size( y );
nlx = prod( szx );
try
    temp = reshape( 1 : nlx, szx );
    ndx_x = builtin( 'subsasgn', temp, S, zeros( szy ) );
catch
    error( lasterr );
end
szx_n = size( ndx_x );

%
% Assign data
%

x = cvx( x );
bx = x.basis_;
if any( szx_n < szx ),
    bx = bx( :, ndx_x );
else
    if any( szx_n > szx ),
        bx( :, end + 1 ) = 0;
        ndx_x( ndx_x == 0 ) = size( bx, 2 );
        bx = bx( :, ndx_x );
        temp = reshape( 1 : prod( szx_n ), szx_n );
    end
    ndx_x = builtin( 'subsref', temp, S );
    ndx_x = ndx_x( : );
    nlz = length( ndx_x );
    y = cvx( y );
    by = y.basis_;
    nx = size( bx, 1 );
    [ ny, my ] = size( by );
    if nx < ny,
        if issparse( by ) & ~issparse( bx ), bx = sparse( bx ); end
        bx( ny, : ) = 0;
    elseif nx > ny,
        if issparse( bx ) & ~issparse( by ), by = sparse( by ); end
        by( nx, : ) = 0;
    end
    if my < nlz,
        by = by( :, ones( 1, nlz ) );
    end
    bx( :, ndx_x ) = by;
end

%
% Create the new object
%

x = cvx( szx_n, bx );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
