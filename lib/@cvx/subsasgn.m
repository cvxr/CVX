function x = subsasgn( x, S, y, cheat )
error( cvx_verify( x, y ) );
if nargin < 4, cheat = []; end
[ prob, x, y ] = cvx_operate( [], x, y );

%
% Check the subscript argument
%

if ~isempty( cheat ),
    err = cvx_subsref_check( nargin, 3, S );
    if ~isempty( err ), error( err ); end
    cheat;
end

%
% Process first (and possibly second) subscript
%

if isnumeric( x ),
    x = cvx( prob, x );
end
switch S(1).type,
    case '()',
        if ~isempty( S(1).subs ),
            szx = size( x );
            szy = size( y );
            nlx = prod( szx );
            if all( szy == 0 ),
                ndx_x = reshape( 1 : nlx, szx );
                try,
                    ndx_x( S(1).subs{:} ) = [];
                catch,
                    error( lasterr );
                end
                szx_n = size( ndx_x );
                bx = cvx_basis( x );
                bx = bx( ndx_x, : );
            else,
                ndx_x = reshape( 1 : nlx, szx );
                try,
                    ndx_x( S(1).subs{:} ) = zeros( szy );
                catch,
                    error( lasterr );
                end
                szx_n = size( ndx_x );
                bx = cvx_basis( x );
                if ~isequal( szx, szx_n ),
                    bx( end + 1, : ) = 0;
                    ndx_x( ndx_x == 0 ) = nlx + 1;
                    bx = bx( ndx_x, : );
                    nlx = prod( szx_n );
                    szx = szx_n;
                end
                ndx_x = reshape( 1 : nlx, szx );
                ndx_x = ndx_x( S(1).subs{:} );
                ndx_x = ndx_x( : );
                by = cvx_basis( y );
                dx = size( bx, 2 );
                dy = size( by, 2 );
                dz = max( dx, dy );
                if dx < dz, bx( end, dz ) = 0; end
                if dy < dz, by( end, dz ) = 0; end
                if size( by, 1 ) < length( ndx_x ),
                    bx( ndx_x, : ) = by( ones( length( ndx_x ), 1 ), : );
                else
                    bx( ndx_x, : ) = by;
                end
            end
            x = cvx( prob, szx_n, bx );
        end
        S(1) = [];
    case '.',
        error( 'All fields are read-only for cvx variable objects.' );
    case '{}',
        error( 'Cell contents reference from a non-cell array object.' );
    otherwise,
        error( 'Invalid subscript structure.' );
end

%
% Process remaining subscripts
%

if ~isempty( S ),
    error( 'Extraneous subscript fields found.' );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
