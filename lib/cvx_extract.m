function [ dbcA, cones, dir, Q, P, exps, dualized ] = cvx_extract( config )

global cvx___
try
    pstr = cvx___.problems(end);
catch
    error( 'CVX:NoModel', 'No CVX model is present.' );
end

%%%%%%%%%%%%%%
% EXTRACTION %
%%%%%%%%%%%%%%

dbcA = pstr.objective;
if pstr.geometric,
    dbcA = log( dbcA );
end
if numel( dbcA ) > 1,
    dbcA = cvx( [1,1], sum( cvx_basis( dbcA ), 2 ) );
end
nold = length( cvx___.classes );
if isempty( dbcA ),
    dir = 1;
    dbcA = cvx( [ 1, 1 ], [] );
elseif strcmp( pstr.direction, 'minimize' ) || strcmp( pstr.direction, 'epigraph' ),
    dbcA = -dbcA;
    dir = 1;
else
    dir = -1;
end
dbcA = cvx_basis( dbcA );
dbcA( end + 1 : nold, : ) = 0;

% Equality constraints
AA = cvx___.equalities;
if pstr.checkpoint(2) > 0,
    npre = sum( cellfun( @(x)size(x,2), AA(1:pstr.checkpoint(2)) ) );
    AA   =  AA( pstr.checkpoint(2) + 1 : end, : );
else
    npre = 0;
end
AA   = cellfun( @(x)vertcat(sparse(x),sparse(nold-size(x,1),size(x,2))), AA, 'UniformOutput', false );
dbcA = horzcat( dbcA, AA{:} );
mold = size(dbcA,2) + npre;
clear AA

% Nonlinearities
used    = full( any( dbcA, 2 ) );
used(1) = true;
rsv     = zeros( nold, 1 ); 
rsv(1)  = 1;
u_int   = false;
nrsv    = 0;
u_exp   = false;
cones   = [];
for k = pstr.checkpoint(3)+1 : length(cvx___.cones),
    cone = cvx___.cones(k);
    ctyp = cone.type;
    ndxs = cone.indices;
    temp = any( reshape( used( ndxs ), size( ndxs ) ), 1 );
    if ~any( temp ), continue; end
    ndxs = ndxs( :, temp );
    used( ndxs ) = true;
    rsv( ndxs ) = 1;
    nnnv = numel(ndxs);
    nrsv = nrsv + nnnv;
    if nargin == 3 && ~isempty( config ),
        switch ctyp,
            case 'exponential',
                if ~config.dualize,
                    error( 'CVX:IncompatibleSolver', 'This solver does not support exponential cones.' );
                end
                nexp = nexp + nnnv / 3;
            case 'semidefinite',
                if ~config.capableSDP && ( size(cones(k).indices,1) > 3 || ~config.capableSOCP ),
                    error( 'CVX:IncompatibleSolver', 'This solver does not support semidefinite cones.' );
                end
            case 'hermitian-semidefinite',
                if ~config.capableSDP && ( size(cones(k).indices,1) > 4 || ~config.capableSOCP ),
                    error( 'CVX:IncompatibleSolver', 'This solver does not support semidefinite cones.' );
                end
            case { 'i_integer', 'i_binary' },
                u_int = true;
                if ~config.INTcapable,
                    error( 'CVX:IncompatibleSolver', 'This solver does not support integer constraints.' );
                end
        end
    else
        if isequal( ctyp(1:2), 'i_' ), u_int = true; end
        if isequal( ctyp, 'exponential' ), u_exp = true; end
    end
    cones = cvx_pushcone( cones, cone.type, ndxs );
end

% Exponentials
nadd = 0;
exps = cvx___.exponential;
if nnz( exps ),
    esrc = find( exps );
    edst = exps( esrc );
    tt   = used(esrc) & used(edst);
    tn   = used(esrc) - used(edst);
    tq   = tn ~= 0;
    exps = [ esrc(tq,:), edst(tq,:), tn(tq,:) ];
    if any( tt ),
        if ~config.dualize,
            error( 'CVX:IncompatibleSolver', 'This solver does not support exponential cones.' );
        end
        % Determine the indices of the exponentials
        u_exp = true;
        esrc  = esrc(tt);
        edst  = edst(tt);
        nexp  = length(esrc);
        nexp3 = 3 * nexp;
        nadd  = nadd + nexp3;
        nrsv  = nrsv + nexp3;
        % Create the exponential cones
        ndim  = reshape( nold+1:nold+nexp3, 3, nexp );
        % Expand Q, P, dbcA
        dbcA(end+nexp3,end) = 0;
        % Add equality consraints to tie the exponential cones to esrc and edst
        % and set the exponential perspective variable to 1
        ndxc = reshape( 1 : 3 * nexp, 3, nexp );
        dbcA = [ dbcA, sparse( ...
            [ esrc(:)' ; ones(1,nexp) ; edst(:)' ; ndim ], ...
            [ ndxc ; ndxc ], ... 
            [ ones(3,nexp) ; -ones(3,nexp) ] ) ];
        cones = cvx_pushcone( cones, 'exponential', ndim );
    end
    if ~isempty( exps ),
        % Look for variables where the exp is used but not the logarithm;
        % these variables must be constrained to be nonnegative.
        tt = exps(:,3) < 0;
        if any( tt ),
            epos = exps(tt,2);
            nlog = length(epos);
            nrsv = nrsv + nlog;
            nadd = nadd + nlog;
            ndim = size(dbcA,1) + (1:nlog);
            dbcA(end+nlog,end) = 0;
            dbcA = [ dbcA, sparse( [epos,ndim']', repmat(1:nlog,[2,1]), ...
                repmat([-1;1],[1,nlog]) ) ];
            cones = cvx_pushcone( cones, 'nonnegative', ndim );
        end
    end
else
    exps = [];
end

if u_exp && u_int,
    error( 'CVX:IncompatibleSolver', 'Exponential variables and integer variables cannot both be used.' );
end

%%%%%%%%%%%%
% PRESOLVE %
%%%%%%%%%%%%

% For the problem
%
%    minimize c' * x + d
%    s.t.     y : A * x + b == 0
%             x \in K
%
% The Lagrangian is
%   
%  L(x,y,z) = c' * x + d - y' * ( A x + b ) - z' * x
%
%                            [ - d  b' 0 ]   [ 1 ]
%            = - [ 1, x' ] * [ - c  A' I ] * [ y ]
%                                            [ z ]
%
% This function provides a smaller [ d, b' ; c, A' ] with no more nonzeros
% that solves an equivalent problem. The original x and y can be recovered
% from the reduced xx and yy by Q*[1;xx] and P*[1;-yy], respectively.

% Add Lagrangian term x'*z
ncones = length(cones);
if config.dualize,
    PP = cell( ncones, 1 );
    nnew = size(dbcA,1);
    for k = 1 : ncones,
        ndxs    = cones(k).indices;
        [nn,nv] = size(ndxs);
        nnnv    = nn * nv;
        SS      = [];
        switch cones(k).type,
            case 'semidefinite',
                SS = 2 * ones(nn,1);
                SS(cumsum([1,round(0.5*(sqrt(8*nn+1)-1)):-1:2])) = 1;
                SS = sparse(1:nn,1:nn,SS);
            case 'hermitian-semidefinite',
                SS = 2 * ones(nn,1);
                SS(cumsum([1,2*sqrt(nn)-1:-2:2])) = 1;
                SS = sparse(1:nn,1:nn,SS);
            case 'exponential',
                persistent SS_exp %#ok
                if isempty( SS_exp ),
                    SS_exp = sparse([0,-1,0;-1,0,0;0,0,exp(-1)]);
                end
                SS = SS_exp;
        end
        PP{k} = sparse(ndxs,1:nnnv,1,nnew,nnnv);
        if ~isempty( SS ),
            PP{k} = PP{k} * cvx_replicate_structure( SS, nv );
        end
    end
    dbcA = [ dbcA, PP{:} ];
    clear PP
end

%
% Initial Q and P matrices
%

[nnew,mnew] = size(dbcA);
if ~all(used),
    ndxs = find(used);
    rsv  = rsv(used,:);
    used(end+1:nnew) = true;
    dbcA = dbcA(used,:);
else
    ndxs = ( 1 : nold )';
end
nred = size(dbcA,1);
Q = sparse( ndxs, 1 : length(ndxs), 1, nold, nred );
P = sparse( [ 1, npre + 2 : mold ], 1 : mold - npre, 1, mold, mnew );
ndxs  = [ ndxs ; ( nold + 1 : nold + nadd )' ];
rsv   = [ rsv  ; ones( nadd, 1 ) ];
ineqs = zeros( 1, mnew );
if config.dualize,
    ineqs(end-nrsv+1:end) = 1;
end
if isempty( dbcA ),
    return;
end

dualized = false;
while true,
    
    last_success = 1;
    while true,
        %
        % STEP 1: Look for columns which differ only by a constant factor.
        % These correspond to redundant equality constraints. These occur
        % often enough as as consequence of our tranformation method, and
        % they cause problems in solvers, so we must eliminate them. Of
        % course, if there are more complex linear dependencies in the
        % equality constraints, we can't do anything about that.
        %
        mcur = size(dbcA,2);
        [ endx, scls ] = cvx_bcompress_mex( dbcA, 0, 1, mcur-nrsv );
        if ~all( diff( endx ) == 1 ),
            last_success = 1;
            mndx = 1 : mcur;
            t2 = endx == mndx;
            xR = sparse( endx, mndx, scls );
            xR = xR(t2,:);
            P = P * cvx_invert_structure( xR );
            dbcA = dbcA(:,t2);
            ineqs = ineqs(t2);
        elseif last_success == 2,
            break;
        end
        %
        % STEP 2: Look for variables that we can eliminate without
        % increasing fill-in. This means looking for rows or columns
        % with only 1, 2, or (in some cases) 3 nonzeros.
        %
        success = false;
        while true,
            [ rows, cols ] = cvx_eliminate_mex( dbcA, 1, rsv, ineqs );
            if ~any( rows ), break; end
            success = true;
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
            A11  = dbcA( rowX, colX );
            A12  = dbcA( rowX, cols );
            A21  = dbcA( rows, colX );
            A22  = dbcA( rows, cols );
            if ( size( A22, 1 ) ~= size( A22, 2 ) || nnz( A22 ) ~= size( A22, 1 ) ),
                error( 'There seems to be an error in the CVX presolver routine.\nPlease report this to the authors; and if possible, include the\ncvx model and data that gave you this error.', 1 ); %#ok
            end
            [ ii, jj, vv ] = find( A22 );
            A22i  = sparse( jj, ii, 1.0 ./ vv );
            temp  = - A22i * A21;
            P     = P( :, colX ) + P( :, cols ) * temp;
            temp  = - A12 * A22i;
            Q     = Q( :, rowX ) + Q( :, rows ) * temp';
            dbcA  = A11 + temp * A21;
            rsv   = rsv( rowX );
            ndxs  = ndxs( rowX );
            ineqs = ineqs( colX );
        end
        if success,
            last_success = 2;
        elseif last_success == 1,
            break;
        end
    end
    %
    % STEP 3: Check to see if dualization will result in a smaller problem.
    % We use the number of equality constraints in each case as our metric.
    % For the primal problem, we have the exact count, but for the dual
    % problem, we may be able to reduce the number by running another pass
    % or two of cvx_eliminate_mex. Here we try running it just once. Thus
    % there is a chance this won't give us the smallest option, but it 
    % should work most of the time. If we want to be very aggressive we
    % can safe off the primal formulation and bring it back in if needed.
    %
    if dualized || ~config.dualize,
        break
    end
    rsv(2:end) = 0;
    [ rows, cols ] = cvx_eliminate_mex( dbcA, 1, rsv, ineqs ); %#ok
    if size(dbcA,2) - nrsv > size(dbcA,1) - nnz(rows),
        dbcA = -dbcA';
        dir  = -dir;
        tmp  = Q; Q = P; P = tmp;
        nnew = size(dbcA,1);
        dualized = true;
        n = nnew - nrsv;
        rsv = zeros(n,1); rsv(1) = 1;
        for k = 1 : length(cones),
            [nn,nv] = size(cones(k).indices);
            nnnv = nn * nv;
            ndxs = n + 1 : n + nnnv;
            rsv(ndxs) = 1; 
            cones(k).indices = reshape(ndxs,nn,nv);
            n = n + nnnv;
        end
        ndxs = 1 : nnew;
    else
        dbcA = dbcA(:,1:end-nrsv);
        P = P(:,1:end-nrsv);
        break;
    end
    
end

%
% STEP 4: Look for redundant slack variables.
%

slacks = false(nnew,1);
for k = 1 : ncones,
    slacks(cones(k).indices(cones(k).slacks>0,:)) = true;
end
slacks = slacks(ndxs);
t_slack = bsxfun( @times, dbcA, slacks );
t_slack = bsxfun( @times, t_slack, ( sum( t_slack ~= 0, 2 ) == 1 ) & ( t_slack(:,1) == 0 ) );
if nnz( t_slack ),
    c_slack = max(0,sum(t_slack>0,1)-1)-max(0,sum(t_slack<0,1)-1);
    c_slack(1) = 0;
    if nnz( c_slack ),
        nneg = zeros(nnew,1);
        if ncones > 1 && isequal( cones(1).type, 'nonnegative' ),
            nneg(cones(1).indices) = 1;
        end
        nneg = nneg(ndxs);
        [ rx, cx ] = find( max( 0, bsxfun( @times, nneg, bsxfun( @times, t_slack, c_slack ) ) ) );
        if ~isempty( rx ),
            dx = diff([0;find(diff(cx))]);
            dx = cumsum(1+sparse(1,cumsum(dx)+1,-dx,1,length(cx)));
            rx = rx(dx<=abs(c_slack(cx)));
            Q(:,rx) = [];
            dbcA(rx,:) = [];
            ndxs(rx) = [];
        end
    end
end
    
%
% Move the cone indices to their new locations and transform
%

nfin = size(dbcA,1);
ndxi = zeros(nnew,1);
ndxi(ndxs) = 1 : length(ndxs);
cone2 = []; Pr = []; Pc = []; Pv = [];
for k = 1 : length(cones),
    cone    = cones(k);
    ctype   = cone.type;
    ndxs    = cone.indices;
    nn      = size(ndxs,1);
    ndxs    = nonzeros(ndxi(ndxs));
    nv      = numel(ndxs) / nn;
    ndxs    = reshape( ndxs, nn, nv );
    switch cone.type,
        case 'lorentz',
            if nn == 2 && config.preferLP,
                beta = sqrt(0.5);
                ndx1 = ndxs(1,:);
                ndx2 = ndxs(2,:);
                Pr = [Pr,ndx1,ndx1,ndx2,ndx2]; %#ok
                Pc = [Pc,ndx1,ndx2,ndx1,ndx2]; %#ok
                Pv = [Pv,beta*[-ones(1,nv),ones(1,3*nv)]]; %#ok
                ctype = 'nonnegative';
            end
        case 'rotated-lorentz',
            if config.unrotateSOCP,
                ndx1 = ndxs(end-1,:);
                ndx2 = ndxs(end,:);
                beta = sqrt(0.5);
                Pr = [Pr,ndx1,ndx1,ndx2,ndx2]; %#ok
                Pc = [Pc,ndx1,ndx2,ndx1,ndx2]; %#ok
                Pv = [Pv,beta*[-ones(1,nv),ones(1,3*nv)]]; %#ok
                ctype = 'lorentz';
            end
        case 'semidefinite',
            if nn == 3 && config.preferSOCP,
                ndx1 = ndxs(1,:);
                ndx2 = ndxs(2,:);
                ndx3 = ndxs(3,:);
                if config.unrotateSOCP,
                    Pr = [Pr,ndx1,ndx1,ndx2,ndx3,ndx3]; %#ok
                    Pc = [Pc,ndx2,ndx3,ndx1,ndx2,ndx3]; %#ok
                    Pv = [Pv,-ones(1,nv),ones(1,4*nv)]; %#ok
                    ctype = 'lorentz';
                else
                    beta = sqrt(0.5);
                    Pr = [Pr,ndx1,ndx2,ndx3]; %#ok
                    Pc = [Pc,ndx2,ndx1,ndx3]; %#ok
                    Pv = [Pv,ones(1,nv)/beta,ones(1,nv),ones(1,nv)/beta]; %#ok
                    ctype = 'rotated-lorentz';
                end
            end
        case 'hermitian-semidefinite',
            if nn == 4 && config.preferSOCP,
                ndx1 = ndxs(1,:);
                ndx2 = ndxs(2,:);
                ndx3 = ndxs(3,:);
                ndx4 = ndxs(4,:);
                if config.unrotateSOCP,
                    Pr = [Pr,ndx1,ndx1,ndx2,ndx3,ndx4,ndx4]; %#ok
                    Pc = [Pc,ndx3,ndx4,ndx1,ndx2,ndx3]; %#ok
                    Pv = [Pv,-ones(1,nv),ones(1,4*nv)]; %#ok
                    ctype = 'lorentz';
                else
                    beta = sqrt(0.5);
                    Pr = [Pr,ndx1,ndx2,ndx3,ndx4]; %#ok
                    Pc = [Pc,ndx3,ndx1,ndx2,ndx4]; %#ok
                    Pv = [Pv,ones(1,nv)/beta,ones(1,2*nv),ones(1,nv)/beta]; %#ok
                    ctype = 'rotated-lorentz';
                end
            elseif config.convertCSDP,
                n2 = sqrt( nn );
                [rs,cs,vs] = find(cvx_create_structure([n2,n2],'hermitian'));
                tt = rem(cs-1,n2) >= floor((cs-1)/n2);
                tr = real(vs) ~= 0;
                is = rs(tt&~tr);
                rs = rs(tt&tr);
                [rd,cd] = find(cvx_create_structure(2*[n2,n2],'symmetric'));
                cc  = floor((cd-1)/(n2*2));
                rr  = rem(cd-1,n2*2);
                rd1 = rd(rr<n2&cc<n2&rr>=cc);
                rd2 = rd(rr>=n2&cc>=n2&rr>=cc);
                id1 = rd(rr>=n2&rr>cc+n2);
                tt  = rr>=n2&cc<n2&rr<cc+n2;
                id2 = rd(tt);
                [dummy,tt] = sort(cc(tt)+n2*rr(tt)); %#ok
                id2 = id2(tt);
                ndxs = [ndxs;reshape(nfin+1:nfin+n2*(n2+1)*nv,[],nv)]; %#ok
                Pr = [ Pr, vec( [ ndxs(rs,:)  ; ndxs(rs,:)  ; ndxs(is,:)  ; ndxs(is,:)  ] )' ]; %#ok
                Pc = [ Pc, vec( [ ndxs(rd1,:) ; ndxs(rd2,:) ; ndxs(id1,:) ; ndxs(id2,:) ] )' ]; %#ok
                Pv = [ Pv, vec( [ ones(2*numel(rs)+numel(is),nv) ; -ones(numel(is),nv) ] )' ]; %#ok
                ctype = 'semidefinite';
                nfin = ndxs(end);
            end
    end
    cone2 = cvx_pushcone( cone2, ctype, ndxs );
end
cones = cone2;
if ~isempty(Pr),
    Pd = true(1,size(dbcA,1),1);
    Pd(Pr) = false;
    Pd = find(Pd);
    Pr = [Pr,Pd];
    Pc = [Pc,Pd];
    Pv = [Pv,ones(size(Pd))];
    Pd = sparse(Pr,Pc,Pv);
    dbcA = Pd' * dbcA;
    Q    = Q * Pd;
end
            
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
