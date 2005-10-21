function eliminate( prob )

global cvx___
prob = index( prob );
p = cvx___.problems( prob );
if p.locked, return; end

rsv = int32( p.reserved( : )' );
nn  = length( p.reserved );
A   = cvx_basis( p.equalities );
if size( A, 2 ) < nn, A( :, nn ) = 0; end
nobj = prod( size( p.objective ) );
preserve_dual = nobj <= 1 & ~isempty( p.duals );
if nobj > 0,
    c = cvx_basis( p.objective );
    if size( c, 2 ) < nn, c( :, nn ) = 0; end
    if isequal( p.direction, 'minimize' ), c = -c; end
    if ~p.complete,
        % For incomplete problems, we want to preserve the exact form
        % of the objective function(s), each which consists of a single
        % variable with a +1 coefficient. This is necessary to insure
        % that convexity verification functions properly. Once the
        % problem has been fully integrated into the parent, these
        % variables can be elimintaed (or not) as needed.
        rsv( any( c, 1 ) ) = 1;
    end
else,
    c = sparse( [], [], [], 1, nn );
    nobj = 1;
end

mm = size( A, 1 );
if preserve_dual,
    P  = speye( nobj + mm );
end
Q  = speye( nn );
ndxs = [ 1 : nn ]';

iter     = 0;
reduced  = 1;
success  = false;

%
% Eliminate variables that are completely unused
%

cols = ~( any( A, 1 ) | any( c, 1 ) | rsv ~= 0 );
if any( cols ),
    warning( 'Unused variables found and eliminated.' );
    ce = find( ~cols );
    A  = A( :, ce );
    c  = c( :, ce );
    Q  = Q( :, ce );
    Q( cols, 1 ) = NaN;
    rsv  = rsv( ce );
    ndxs = ndxs( ce );
    success = true;
end

while 1,
    iter = iter + 1;
    if reduced == iter - 2,
        break;
    end
    if rem( iter, 2 ) == 1,
        [ xL, A, xI ] = cvx_bcompress( A );
        if size( xL, 1 ) == size( xL, 2 ),
            continue;
        end
        if preserve_dual,
            P = [ P( :, 1 : nobj ), P( :, nobj + 1 : end ) * xI ];
        end
        Az = A ~= 0;
        reduced = iter;
        success = true;
    else,
        while ~isempty( A ),
            %
            % Select columns to eliminate. More about this later. In short,
            % we're looking for an invertible diagonal block of A whose
            % rows and columns are sparse, so that no fill-in occurs.
            %

            [ rows, cols, rsv ] = cvx_eliminate_mex( A, c, rsv );
            if ~any( rows ), break; end
            rows = rows ~= 0;
            cols = cols ~= 0;
            rowX = ~rows;
            colX = ~cols;

            %
            % [ c1  c2  ] = [ I 0   0 ] [ c1-c2*A22i*A21   0    ] [ I   0   ]
            % [ A11 A12 ]   [ 0 I A12 ] [ A11-A12*A22i*A21 0    ] [ A21 A22 ]
            % [ A21 A22 ] = [ 0 0 A22 ] [ 0                A22i ]
            %

            c1   = c( :,    colX );
            c2   = c( :,    cols );
            A11  = A( rowX, colX );
            A12  = A( rowX, cols );
            A21  = A( rows, colX );
            A22  = A( rows, cols );
            [ ii, jj, vv ] = find( A22 );
            A22i = sparse( jj, ii, 1.0 ./ vv );
            if preserve_dual,
                temp = - A22i' * [ c2 ; A12 ]';
                tndx = find(rowX) + nobj;
                P    = P( :, [ 1 : nobj, tndx(:)' ] ) + P( :, find(rows) + nobj ) * temp;
            end
            temp = - A22i * A21;
            Q    = Q( :, colX ) + Q( :, cols ) * temp;
            A    = A11 + A12 * temp;
            c    = c1  + c2  * temp;
            rsv  =  rsv( colX );
            ndxs = ndxs( colX );
            reduced = iter;
            success = true;
        end
    end
end

if ~success,
    return
end

if ~isempty( p.objective ),
    if isequal( p.direction, 'minimize' ), c = -c; end
    p.objective = cvx( p.self, size( p.objective ), c );
end
p.equalities = cvx( p.self, size( A, 1 ), A );
p.variables = apply_map( p.variables, Q );
if preserve_dual,
    p.duals = apply_map( p.duals, P );
else,
    p.duals = [];
end
ndxs2 = zeros( 1, length( p.reserved ) );
ndxs2( ndxs ) = 1;
ndxs2 = cumsum( ndxs2 );
for k = 1 : length( p.cones ),
    temp = p.cones( k ).indices;
    p.cones( k ).indices = reshape( ndxs2( temp ), size( temp ) );
end
p.reserved = p.reserved( ndxs );
p.vexity   = p.vexity( ndxs );
p.locked   = true;

cvx___.problems( prob ) = p;

function obj = apply_map( obj, map )

s = size( obj );
n = prod( s );
switch class( obj ),
    case 'cell',
        for k = 1 : n,
            obj{ k } = apply_map( obj{ k }, map );
        end
    case 'struct',
        obj = cell2struct( apply_map( struct2cell( obj ), map ), fieldnames( obj ), 1 );
    otherwise,
        bobj = cvx_basis( obj );
        obj  = cvx( problem( obj ), s, bobj * map( 1 : size( bobj, 2 ), : ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
