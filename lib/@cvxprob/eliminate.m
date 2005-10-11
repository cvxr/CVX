function eliminate( prob )

global cvx___
prob = index( prob );
p = cvx___.problems( prob );
if p.locked, return; end

rsv = p.reserved( : )';
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
        rsv( any( c, 1 ) ) = true;
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
b_sparse = true;
c_sparse = true;
reduced  = 1;
success  = false;

%
% Eliminate variables that are completely unused
%

cols = ~( any( A, 1 ) | any( c, 1 ) | rsv );
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
        while 1,
            %
            % Select columns to eliminate. More about this later. In short,
            % we're looking for an invertible diagonal block of A whose
            % rows and columns are sparse, so that no fill-in occurs.
            %
            [ mm, nn ] = size( A );
            if mm == 0, break; end
            [ rowA, colA ] = find( A );
            [ rowC, colC ] = find( c );
            Ac = full( sparse( colA, 1, 1, nn, 1 ) );
            if c_sparse, Ac = Ac + full( sparse( colC, 1, 1, nn, 1 ) ); end
            Ac( rsv ) = Inf;
            Ar = full( sparse( rowA, 1, 1, mm, 1 ) );
            if ~b_sparse, Ar = Ar - full( A( :, 1 ) ~= 0 ); end
            rmax = max( Ar ) + 1;
            rcnt = sparse( rowA, colA, Ar(rowA,:)-rmax*~rsv(:,colA)', mm, nn );
            [ rcnt, rows ] = min( rcnt, [], 1 );
            rcnt = rcnt' + rmax;
            cols = find( Ac <= 2 | rcnt <= 2 | ( rcnt + Ac <= 7 ) );
            rows = rows( cols );
            %
            % Insure that the submatrix we have selected is diagonal. If it 
            % is not, then it means that a variable that is being
            % eliminated is serving as another's replacement. That cycle
            % must be eliminated before we can proceed.
            %
            A22 = A( rows, cols );
            if nnz( A22 ) > size( A22, 1 ),
                [ rd, cd ] = find( A22 );
                temp = rd == cd;
                rd( temp ) = []; cd( temp ) = [];
                temp = max( [ rd( : ), cd( : ) ], [], 2 );
                rows( temp ) = []; 
                cols( temp ) = [];
                A22 = A( rows, cols );
                if nnz( A22 ) > size( A22, 1 ),
                    temp = sum( A22 ~= 0, 1 ) == 1;
                    rows = rows( temp ); 
                    cols = cols( temp );
                    A22 = A( rows, cols );
                end
            end
            if isempty( cols ),
                break;
            end
            %
            % [ c1  c2  ] = [ I 0   0 ] [ c1-c2*A22i*A21   0 ]
            % [ A11 A12 ]   [ 0 I A12 ] [ A11-A12*A22i*A21 0 ]
            % [ A21 A22 ] = [ 0 0 A22 ] [ A22i*A21         I ]
            %
            nr   = length( cols );
            rowX = 1 : mm; rowX( rows ) = [];
            colX = 1 : nn; colX( cols ) = [];
            c1   = c( :,    colX );
            c2   = c( :,    cols );
            A11  = A( rowX, colX );
            A12  = A( rowX, cols );
            A21  = A( rows, colX );
            A22i = 1 : nr;
            A22i = sparse( A22i, A22i, 1.0 ./ diag( A22 ), nr, nr );
            if preserve_dual,
                temp = - A22i * [ c2 ; A12 ]';
                P    = P( :, [ 1 : nobj, rowX + nobj ] ) + P( :, rows + nobj ) * temp;
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
