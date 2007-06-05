function [ dbcA, cones, objsize, dir, Q, P, esrc, edst ] = extract( p, destructive )
if nargin < 2, destructive = false; end

global cvx___
p = cvx___.problems( index( p ) );

%
% Sizes
%

nb = length( cvx___.reserved );
nres = length( p.t_variable );
if nres == 1 | all( p.t_variable ),
    need_used = false;
    n = nb;
else
    need_used = true;
    used = p.t_variable;
    used( end + 1 : nb, : ) = true;
    offset = cumsum( used );
    n = offset( end );
end

%
% Objective
%

dbcA    = p.objective;
objsize = size( dbcA );
nobj    = prod( objsize );
if nobj == 0,
    nobj = 1;
    dir = 1;
    objsize = [ 1, 1 ];
    dbcA = cvx( [ 1, 1 ], [] );
elseif strcmp( p.direction, 'minimize' ) | strcmp( p.direction, 'epigraph' ),
    dbcA = vec( dbcA );
    dir = 1;
else
    dbcA = -vec( dbcA );
    dir = -1;
end

%
% Equality constraints
%

AA = cvx___.equalities;
ineqs = cvx___.needslack;
if p.n_equality > 0,
    AA = AA( p.n_equality + 1 : end, : );
    ineqs = ineqs( p.n_equality + 1 : end, : );
    if destructive,
        cvx___.equalities( p.n_equality + 1 : end ) = [];
        cvx___.needslack( p.n_equality + 1 : end ) = [];
    end
elseif destructive,
    cvx___.equalities = cvx( [ 0, 1 ], [] );
    cvx___.needslack = logical( zeros( 0, 1 ) );
end
if ~isempty( AA ),
    ineqs = [ logical(zeros(nobj,1)) ; ineqs ];
    dbcA = [ dbcA ; AA ];
    clear AA
end

%
% Linear forms
%

if p.n_linform > 0,
    A1 = cvx___.linforms( p.n_linform + 1 : end, : );
    A2 = cvx___.linrepls( p.n_linform + 1 : end, : );
    if destructive,
        cvx___.linforms( p.n_linform + 1 : end ) = [];
        cvx___.linrepls( p.n_linform + 1 : end ) = [];
    end
else
    A1 = cvx___.linforms;
    A2 = cvx___.linrepls;
    if destructive,
        cvx___.linforms = cvx( [ 0, 1 ], [] );
        cvx___.linrepls = cvx( [ 0, 1 ], [] );
    end
end
if ~isempty( A1 ),
    zV = cvx_vexity( A2 ); 
    zQ = ( zV == 0 ) - zV;
    dbcA = [ dbcA ; minus( zQ .* A1, zQ .* A2, true ) ];
    ineqs( end + 1 : end + length( A1 ), : ) = zV ~= 0;
    clear A1 A2 zV zQ
end

%
% Univariable forms
%

if p.n_uniform > 0,
    A1 = cvx___.uniforms( p.n_uniform + 1 : end, : );
    A2 = cvx___.unirepls( p.n_uniform + 1 : end, : );
    if destructive,
        cvx___.uniforms( p.n_uniform + 1 : end ) = [];
        cvx___.unirepls( p.n_uniform + 1 : end ) = [];
    end
else
    A1 = cvx___.uniforms;
    A2 = cvx___.unirepls;
    if destructive,
        cvx___.uniforms = cvx( [ 0, 1 ], [] );
        cvx___.unirepls = cvx( [ 0, 1 ], [] );
    end
end
if ~isempty( A2 ),
    zV = cvx_vexity( A2 );
    zQ = ( zV == 0 ) - zV;
    dbcA = [ dbcA ; minus( zQ .* A1, zQ .* A2, true ) ];
    ineqs( end + 1 : end + length( A1 ), : ) = zV ~= 0;
    clear A1 A2 zV zQ
end

%
% Convert to basis
%

dbcA = cvx_basis( dbcA );
nA = size( dbcA, 1 );
if nA < nb,
    dbcA( end + 1 : nb, : ) = 0;
elseif nb < nA,
    dbcA = dbcA( 1 : nb, : );
end

%
% Convert inequalities to equalities if they have available slack
%

if any( ineqs ),
    slacks = cvx___.canslack;
    if any( slacks ),
        ndxs   = find( ineqs );
        sterms = dbcA( slacks, ineqs );
        oterms = dbcA( slacks, 1 : nobj );
        if nnz( sterms ),
            sdirec = cvx___.vexity( slacks );
            pslack = sum( sterms < 0, 2 ) == 1 & sdirec >= 0 & ~any( oterms > 0, 2 );
            nslack = sum( sterms > 0, 2 ) == 1 & sdirec <= 0 & ~any( oterms < 0, 2 );
            qslack = pslack & nslack;
            temp = any( sterms( pslack & ~nslack, : ) < 0, 1 ) | ...
                   any( sterms( nslack & ~pslack, : ) > 0, 1 );
            ineqs( ndxs( temp ) ) = false;
            if any( qslack ),
                ndxs = ndxs( ~temp );
                [ rr, cc, vv ] = find( sterms( qslack, ~temp ) );
                [ c1, ci ] = unique( cc ); [ c2, ri ] = unique( rr( ci ) );
                [ c2, rj ] = unique( rr ); [ c2, cj ] = unique( cc( rj ) );
                if length( c2 ) < length( ri ), c2 = c1( ri ); end
                ineqs( ndxs( c2 ) ) = false;
            end
        end
    end
end

%
% Eliminate unused variables
%

if need_used,
    dbcA = dbcA( used, : );
end

%
% Add slack variables, part 1
%

nsl = nnz( ineqs );
if nsl ~= 0,
    sndxs = size( dbcA, 1 ) + [ 1 : nsl ];
    dbcA = [ dbcA ; sparse( 1 : nsl, find( ineqs ), -1, nsl, length( ineqs ) ) ];
end

%
% Cones
%

if nargout < 3 & ~destructive,
    return
end
if need_used,
    cones = [];
    ocones = [];
    for k = 1 : length( cvx___.cones ),
        cone = cvx___.cones( k );
        temp = reshape( used( cone.indices ), size( cone.indices ) );
        temp = any( temp, 1 );
        if any( temp ),
            ncone = cone;
            ncone.indices = ncone.indices( :, temp );
            if need_used,
                ncone.indices = reshape( offset( ncone.indices ), size( ncone.indices ) );
            end
            if isempty( cones ),
                cones = ncone;
            elseif isequal( ncone.type, 'nonnegative' ),
                cones = [ ncone, cones ];
            else
                cones = [ cones, ncone ];
            end
            if destructive,
                temp = ~temp;
                if any( temp ),
                    ncone = cone;
                    ncone.indices = ncone.indices( :, temp );
                    if isempty( ocones ),
                        ocones = ncone;
                    else
                        ocones = [ ocones, ncone ];
                    end
                end
            end
        end
    end
    if destructive,
        cvx___.cones = ocones;
    end
else
    cones = cvx___.cones;
    if destructive,
        cvx___.cones = [];
    end
end

%
% Add slack variables, part 2
%

if nsl ~= 0,
    ncone = struct( 'type', 'nonnegative', 'indices', sndxs );
    if isempty( cones ),
        cones = ncone;
    elseif isequal( cones( 1 ).type, 'nonnegative' ),
        cones( 1 ).indices = [ cones( 1 ).indices, ncone.indices ];
    else
        cones = [ ncone, cones ];
    end
end

%
% Exponential and logarithm indices
%

if nargout > 6 | destructive,
    esrc = find( cvx___.exponential );
    edst = full( cvx___.exponential( esrc ) );
    if need_used,
        tt = ~( used( esrc ) & used( edst ) );
        esrc( tt ) = [];
        edst( tt ) = [];
    end
end

%
% Reserved flags
%

if destructive,
    cvx___.reserved( nres + 1 : end ) = [];
    cvx___.geometric( nres + 1 : end ) = [];
    cvx___.exponential( nres + 1 : end ) = [];
    cvx___.logarithm( nres + 1 : end )   = [];
    cvx___.exponential( cvx___.exponential > nres ) = 0;
    cvx___.logarithm( cvx___.logarithm > nres ) = 0;
end

%
% Q and P matrices
%

if nargout > 4,
    if need_used,
        Q = sparse( find( used ), 1 : n, 1, nb, n + nsl );
    else
        Q = speye( nb, nb + nsl );
    end
end

if nargout > 5,
    npre = p.n_equality;
    ntot = length( cvx___.equalities );
    P    = sparse( npre + 1 : ntot, nobj + 1 : nobj + ( ntot - npre ), ...
                   1, ntot, size( dbcA, 2 ) );
end
