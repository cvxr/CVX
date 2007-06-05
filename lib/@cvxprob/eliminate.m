function [ dbCA, cones, objsize, dir, Q, P, esrc, edst ] = eliminate( prob, destructive )
if nargin < 2, destructive = false; end

global cvx___
[ dbCA, cones, objsize, dir, Q, P, esrc, edst ] = extract( prob, destructive );

if ~issparse( dbCA ),
    dbCA = sparse( dbCA );
end
nn = size( dbCA, 1 );
rsv = sparse( 1, 1, 1, nn, 1 );
for k = 1 : length( cones ),
    rsv = rsv + sparse( cones( k ).indices, 1, 1, nn, 1 );
end
if ~isempty( esrc ),
    tt = any( dbCA( esrc, : ), 2 ) & any( dbCA( edst, : ), 2 );
    if any( tt ),
        rsv( esrc( tt ) ) = rsv( esrc( tt ) ) + 1;
        rsv( edst( tt ) ) = rsv( edst( tt ) ) + 1;
    end
end
nobj  = prod( objsize );
rsv   = full( rsv );
ndxs  = [ 1 : nn ]';
nold  = nn;
ondxs = ndxs;

if size( dbCA, 1 ) == nobj,
    return
end

%
% Check for variables that are either completely unused or are found only
% in the objective (and therefore are unbounded). The unused variables are
% eliminated completely. All but *one* of the variables found in the
% objective is eliminated. One is kept around to make sure that the
% unbounded behavior of the problem is preserved.
%

iter = 0;
last_success = 1;
while true,
    
    iter = iter + 1;
    if last_success + 3 == iter, break; end
    if isempty( dbCA ), break; end
    switch rem( iter, 3 ),
        case 1,
            %
            % Look for rows which are completely empty---which implies that the
            % variables are unused---or rows which occur only in a single
            % objective (and no constraints)---which are trivially unbounded.
            % Unbounded variables cannot *all* be eliminated; one must be kept
            % in order to preserve the unbounded behavior. 
            %
            rcnt = sum( dbCA ~= 0, 2 );
            drsv = ~double( rsv );
            for obj = 1 : nobj,
                rows = ( rcnt == ( dbCA( :, obj ) ~= 0 ) ) & drsv;
                if any( rows ),
                    last_success = iter;
                    celm = dbCA( rows, obj );
                    rowX = ~rows;
                    if nnz( celm ),
                        cnrm = norm( celm );
                        keep = Q( :, rows ) * ( celm / cnrm );
                        ndxq = find( rows );
                        ndxq = ndxq( 1 );
                        rowX( ndxq ) = 1;
                    end
                    dbCA = dbCA( rowX, : );
                    rsv  = rsv ( rowX, : );
                    ndxs = ndxs( rowX, : );
                    Q    = Q( :, rowX );
                    if nnz( celm ),
                        dbCA( ndxq, obj ) = cnrm;
                        Q( :, ndxq ) = keep;
                    end
                end
            end
        case 2,
            %
            % Look for columns which differ only by a constant factor; these
            % correspond to redundant equality constraints. These occur often
            % enough, as a consequence of our transformation method, that they
            % need to be identified and eliminated.
            %
            [ xR, dbCA ] = cvx_bcompress( dbCA, 'full', nobj );
            if size( xR, 1 ) ~= size( xR, 2 ), 
                last_success = iter;
                P = P * cvx_invert_structure( xR );
            end
        otherwise,
            %
            % Look for variables that we can eliminate without increasing
            % fill-in. This basically means looking for rows or columns with
            % only 1, 2, or (in some cases) 3 nonzero entries. If we
            % succeed in finding variables, we back the iterator count by
            % one so that we can try again until we fail.
            %
            [ rows, cols ] = cvx_eliminate_mex( dbCA, nobj, rsv );
            if any( rows ),
                last_success = iter;
                iter = iter - 1;
                rows = rows ~= 0;
                cols = cols ~= 0;
                rowX = ~rows;
                colX = ~cols;
                %
                % [ x1^T x2^T ] [ C1 A11 A12 ] [ 1  ]
                %               [ C2 A21 A22 ] [ y1 ] = 0
                %                              [ y2 ]
                %
                % [ x1^T x2^T ] = x1^T [ I -A12*A22i ]
                %
                % [ G Y1^T Y2^T ] = [ G Y1^T ] [ I  0  -C2'*A22i'  ]
                %                              [ 0  I  -A21'*A22i' ]
                %
                A11  = dbCA( rowX, colX );
                A12  = dbCA( rowX, cols );
                A21  = dbCA( rows, colX );
                A22  = dbCA( rows, cols );
                if ( size( A22, 1 ) ~= size( A22, 2 ) | nnz( A22 ) ~= size( A22, 1 ) ),
                    error( sprintf( 'There seems to be an error in the CVX presolver routine.\nPlease report this to the authors; and if possible, include the\ncvx model and data that gave you this error.' ) );
                end
                [ ii, jj, vv ] = find( A22 );
                A22i = sparse( jj, ii, 1.0 ./ vv );
                temp = - A22i * A21;
                P    = P( :, colX ) + P( :, cols ) * temp;
                temp = - A12 * A22i;
                Q    = Q( :, rowX ) + Q( :, rows ) * temp';
                dbCA = A11 + temp * A21;
                rsv  =  rsv( rowX );
                ndxs = ndxs( rowX );
            end
    end        
end

%
% The first nobj columns of P need to be negated in order to produce the right
% numerical results. This is basically because the dcbA matrix should have
% negated the corresponding rows, but we didn't want to. (If you express the
% Lagrangian in terms of the block dcbA matrix this becomes clear.)
%

if ~isempty( P ) & nobj > 0,
    P( :, 1 : nobj ) = - P( :, 1 : nobj );
end

%
% Move the cone indices to their new locations
%

ndxs = full( sparse( ndxs, 1, 1 : length( ndxs ), nold, 1 ) );
for k = 1 : length( cones ),
    cones(k).indices = reshape( ndxs(cones(k).indices), size(cones(k).indices) );
end
if ~isempty( esrc ),
    esrc = ndxs( esrc );
    edst = ndxs( edst );
    tt = esrc & edst;
    esrc = esrc( tt );
    edst = edst( tt );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
