function [ dbcA, cones, dir, Q, P, exps, ineqs ] = extract( pp, destructive )
if nargin < 3 || nargout < 6, doineqs = true; end
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

AA    = cvx___.equalities;
ntot  = cvx___.n_equality;
ineqs = cvx___.needslack;
szs   = cellfun( @(x)size(x,2), AA );
if p.n_equality > 0,
    npre  = sum( szs( 1: p.n_equality ) );
    szs   =    szs( p.n_equality + 1 : end );
    AA    =    AA( p.n_equality + 1 : end, : );
    ineqs = ineqs( p.n_equality + 1 : end, : );
else
    npre  = 0;
end
AA   = cellfun( @(x)vertcat(x,sparse(n-size(x,1),size(x,2))), AA, 'UniformOutput', false );
dbcA = horzcat( dbcA, AA{:} );
ndxs = zeros(1,ntot);
ndxs(cumsum([1,szs(1:end-1)])) = 1;
ineqs = [false;ineqs(cumsum(ndxs),:)];
clear ndxs AA

%
% Determine which inequalities need slack variables. Not all of them do,
% because some may contain variables which themselves can absorb any slack,
% and thus can be converted to equations without sacrificing equivalence.
% We are somewhat conservative in our determinations here: the variable
% must appear *only* in this inequality, and has been identified as free
% to grow without further constraint in the direction of slack. For
% example, consider the inequality
%   x + y <= z
% where y is an epigraph variable of a convex constraint f(w) <= y. If y
% does not *also* appear in the objective coerced against growth, then we
% are free to replace the inequality with the equation
%   x + y == z
% and equivalence is preserved. We attempt to be very conservative here,
% so it is possible that we do not catch all of the cases where
% inequalities may be converted to equations.
%

if any( ineqs ) && any( cvx___.canslack ),
    slacks = find( cvx___.canslack );
    sterms = dbcA( slacks, : );
    oterms = sterms( :, 1 );
    ecount = sum( sterms( :, ~ineqs ) ~= 0, 2 ) - ( oterms ~= 0 );
    sterms = sterms( :, ineqs );
    icount = sum( sterms ~= 0, 2 );
    sterms = sum( sterms, 2 );
    sdirec = vex(cvx___.classes(slacks));
    nslack = icount == 1 & ecount == 0 & sterms .* sdirec <= 0 & sterms .* oterms >= 0;
    ineqs( any( dbcA( slacks( nslack ), : ), 1 ) ) = false;
end

%
% Select the cones used
%

used = full( any( dbcA, 2 ) );
if all( used ),
    cones = cvx___.cones;
else
    cones = [];
    for k = 1 : length( cvx___.cones ),
        cone = cvx___.cones( k );
        temp = any( reshape( used( cone.indices ), size( cone.indices ) ), 1 );
        if any( temp ),
            ncone = cone;
            ncone.indices = ncone.indices( :, temp );
            if isempty( cones ),
                cones = ncone;
            else
                cones = [ cones, ncone ];
            end
        end
    end
end

%
% Add the slack variables
%

if doineqs,
    nsl = nnz( ineqs );
    if nsl ~= 0,
        dbcA = [ dbcA ; sparse( 1 : nsl, find( ineqs ), -1, nsl, length( ineqs ) ) ];
        ncone = struct( 'type', 'nonnegative', 'indices', n+1:n+nsl );
        if isempty( cones ),
            cones = ncone;
        else
            tt = find(strcmp({cones.type},'nonnnegative'));
            if ~isempty( tt ),
                cones(tt(1)).indices = [ cones(tt(1)).indices, ncone.indices ];
            else
                cones = [ ncone, cones ];
            end
        end
    end
    ineqs = [];
else
    ineqs = find(ineqs);
end

%
% Q and P matrices
%

used = find( used );
Q = sparse( used, used, 1, n, n + nsl );
P = sparse( [ 1, npre + 2 : ntot + 1 ], 1 : ntot - npre + 1, 1, ntot + 1, size( dbcA, 2 ) );

%
% Exponential and logarithm indices
%

exps = cvx___.exponential;
if ~isempty( exps ),
    esrc = find( exps );
    edst = exps( esrc );
    tt   = any(dbcA(esrc,:),2) & any(dbcA(edst,:),2);
    tn   = ~tt;
    exps = [ esrc(tn), edst(tn) ];
    if any( tt ),
        % Determine the indices of the exponentials
        esrc = esrc(tt);
        edst = edst(tt);
        nexp = length(esrc);
        nexp3 = 3 * nexp;
        % Create the exponential cones
        ncone.type = 'exponential';
        ncone.indices = reshape( n+nsl+(1:nexp3), 3, nexp );
        % Expand Q, P, dbCA
        Q(end,end+nexp3) = 0;
        P(end,end+nexp3) = 0;
        dbcA(end+nexp3,end) = 0;
        % Add equality consraints to tie the exponential cones to esrc and edst
        % and set the exponential perspective variable to 1
        ndxc = reshape( 1 : 3 * nexp, 3, nexp );
        dbcA = [ dbcA, sparse( ...
            [ esrc(:)' ; ones(1,nexp) ; edst(:)' ; ncone.indices ], ...
            [ ndxc ; ndxc ], ... 
            [ ones(3,nexp) ; -ones(3,nexp) ] ) ];
        if isempty( cones ),
            cones = ncone;
        else
            tt = find(strcmp({cones.type},'exponential'));
            if ~isempty( tt ),
                cones(tt(1)).indices = [ cones(tt(1)).indices, ncone.indices ];
            else
                cones = [ cones, ncone ];
            end
        end
    end
end

%
% Reserved flags
%

if destructive,
    erase( pp );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
