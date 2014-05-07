function [ dbcA, cones, dir, Q, P, exps ] = extract( pp, destructive )
if nargin < 2 || nargout < 5, destructive = false; end

global cvx___
[ pn, p ] = verify( pp ); %#ok
persistent vex
if isempty( vex ),
    vex = [0,0,0,NaN,-1,-1,-1,0,0,0,1,1,1,NaN,NaN,-1,-1,NaN,NaN,-1,NaN,NaN,NaN]';
end
n = length( cvx___.classes );

%
% Objective
%

dbcA = p.objective;
if numel( dbcA ) > 1,
    dbcA = cvx( [1,1], sum( cvx_basis( dbCA ), 2 ) );
end
if isempty(p.objective),
    dir = 1;
    dbcA = cvx( [ 1, 1 ], [] );
elseif strcmp( p.direction, 'minimize' ) || strcmp( p.direction, 'epigraph' ),
    dir = 1;
else
    dbcA = -dbcA;
    dir = -1;
end
dbcA = cvx_basis( dbcA );
dbcA( end + 1 : n, : ) = 0;

%
% Equality constraints
%

AA     = cvx___.equalities;
ntot   = cvx___.n_equality;
szs    = cellfun( @(x)size(x,2), AA );
if p.n_equality > 0,
    npre  = sum( szs( 1: p.n_equality ) );
    AA    =  AA( p.n_equality + 1 : end, : );
else
    npre  = 0;
end
AA   = cellfun( @(x)vertcat(sparse(x),sparse(n-size(x,1),size(x,2))), AA, 'UniformOutput', false );
dbcA = horzcat( dbcA, AA{:} );
clear AA

%
% Nonlinearities
%

used  = full( any( dbcA, 2 ) );
used(1) = true;
cones = cvx___.cones;
ncone = length(cones);
allu  = all(used);
if ~allu,
    prune = false;
    for k = 1 : ncone,
        cone = cones(k);
        ndxs = cone.indices;
        temp = any( reshape( used( ndxs ), size( ndxs ) ), 1 );
        if ~all( temp ),
            ndxs = ndxs( :, temp );
            used( ndxs ) = true;
            if isempty( ndxs ), prune = true; end
            cones( k ).indices = temp;
        end
    end
    if prune,
        cones = cones(cellfun(@(x)~isempty(x.indices),cones));
    end
    allu = all( used );
end

%
% Exponentials
%

exps = cvx___.exponential;
if ~isempty( exps ),
    esrc = find( exps );
    edst = exps( esrc );
    tt   = used(esrc) & used(edst);
    tn   = used(esrc) - used(edst);
    tq   = tn ~= 0;
    exps = [ esrc(tq), edst(tq), tn(tq) ];
    if any( tt ),
        % Determine the indices of the exponentials
        esrc  = esrc(tt);
        edst  = edst(tt);
        nexp  = length(esrc);
        nexp3 = 3 * nexp;
        % Create the exponential cones
        indices = reshape( n+1:n+nexp3, 3, nexp );
        % Expand Q, P, dbCA
        dbcA(end+nexp3,end) = 0;
        % Add equality consraints to tie the exponential cones to esrc and edst
        % and set the exponential perspective variable to 1
        ndxc = reshape( 1 : 3 * nexp, 3, nexp );
        dbcA = [ dbcA, sparse( ...
            [ esrc(:)' ; ones(1,nexp) ; edst(:)' ; indices ], ...
            [ ndxc ; ndxc ], ... 
            [ ones(3,nexp) ; -ones(3,nexp) ] ) ];
        done = false;
        if ~isempty( cones )
            tt = find( strcmp( { cones.type }, 'exponential' ) );
            if ~isempty( tt ),
                cones(tt(1)).indices = [ cones(tt(1)).indices, indices ];
                done = true;
            end
        end
        if ~done,
            cones = [ cones, struct( 'type', 'exponential', 'indices', indices ) ];
        end
    end
end

%
% Q and P matrices
%

if allu,
    Q = sparse( 1 : n, 1 : n, 1, n, size(dbcA,1) );
else
    nq = nnz(used);
    nn = size(dbcA,1);
    Q  = sparse( find(used), 1 : nq, 1, n, nq+(nn-n) );
    used(end+1:nn) = true;
    dbcA   = dbcA(used,:);
    if ~isempty( cones ),
        ndxi = zeros(1,nn);
        used(n+1:nn) = true;
        ndxi(used) = 1:nq+(nn-n);
        for k = 1 : length(cones),
            cones(k).indices = ndxi(cones(k).indices);
        end
    end
end
P = sparse( [ 1, npre + 2 : ntot + 1 ], 1 : ntot - npre + 1, 1, ntot + 1, size(dbcA,2) );

%
% Reserved flags
%

if destructive,
    erase( pp );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
