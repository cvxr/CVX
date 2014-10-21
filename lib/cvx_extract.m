function [ dbcA, cones, dir, Q, P, exps, dualized ] = cvx_extract( config, name )

global cvx___
try
    pstr = cvx___.problems(end);
catch
    cvx_throw( 'No CVX model is present.' );
end

%%%%%%%%%%%%%%
% EXTRACTION %
%%%%%%%%%%%%%%

dualize = isempty( config ) || ( isfield( config, 'dualize' ) && config.dualize );
nold = length( cvx___.classes );
dbcA = pstr.objective;
if isempty( dbcA ),
    dbcA = sparse( nold, 1 );
    dir = 1;
else
    if abs( pstr.direction ) > 1,
        dbcA = log( dbcA );
    end
    if numel( dbcA ) > 1,
        dbcA = cvx( [1,1], sum( cvx_basis( dbcA ), 2 ) );
    end
    if isempty( dbcA ),
        dir = 1;
        dbcA = cvx( [ 1, 1 ], [] );
    elseif pstr.direction > 0,
        dbcA = -dbcA;
        dir = 1;
    else
        dir = -1;
    end
    dbcA = cvx_basis( dbcA );
    dbcA( end + 1 : nold, : ) = 0;
end

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

% Nonlinearities: pass 1
used    = full( any( dbcA, 2 ) );
used(1) = true;
rsv     = zeros( nold, 1 ); 
rsv(1)  = 1;
nrsv    = 0;
u_int   = false;
u_exp   = false;
u_nneg  = false;
cones   = [];
for k = pstr.checkpoint(3)+1 : length(cvx___.cones),
    cone = cvx___.cones(k);
    ndxs = cone.indices;
    temp = any( reshape( used( ndxs ), size( ndxs ) ), 1 );
    if ~any( temp ), continue; end
    switch cone.type,
        case 'nonnegative', u_nneg = true;
        case { 'integer', 'binary' }, u_int = true;
    end
    ndxs = ndxs( :, temp );
    nrsv = nrsv + numel( ndxs );
    used( ndxs ) = true;
    rsv( ndxs ) = 1;
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
        % Determine the indices of the exponentials
        esrc  = esrc(tt);
        edst  = edst(tt);
        nexp  = length(esrc);
        nexp3 = 3 * nexp;
        nadd  = nadd + nexp3;
        nrsv  = nrsv + nexp3;
        u_exp = true;
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
            nadd = nadd + nlog;
            nrsv = nrsv + nlog;
            ndim = size(dbcA,1) + (1:nlog);
            u_nneg = true;
            dbcA(end+nlog,end) = 0;
            dbcA = [ dbcA, sparse( [epos,ndim']', repmat(1:nlog,[2,1]), ...
                repmat([-1;1],[1,nlog]) ) ];
            cones = cvx_pushcone( cones, 'nonnegative', ndim );
        end
    end
else
    exps = [];
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
if dualize,
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
            case 'hermitian_semidefinite',
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
ineqs(1) = 1;
if dualize,
    ineqs(end-nrsv+1:end) = 1;
end
if isempty( dbcA ),
    return;
end

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
    mndx = 1 : mcur;
    t2 = endx == mndx;
    xR = sparse( endx, mndx, scls );
    t2 = t2 & scls;
    xR = xR(t2,:);
    P = P * cvx_invert_structure( xR )';
    dbcA = dbcA(:,t2);
    ineqs = ineqs(t2);
end

dualized = false;
while true,
    %
    % STEP 2: Look for variables that we can eliminate without
    % increasing fill-in. This means looking for rows or columns
    % with only 1, 2, or (in some cases) 3 nonzeros.
    %
    while true,
        [ rows, cols ] = cvx_eliminate_mex( dbcA, 1, rsv, ineqs );
        if ~any( rows ), break; end
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
            cvx_throw( 'There seems to be an error in the CVX presolver routine.\nPlease report this to the authors; and if possible, include the\ncvx model and data that gave you this error.', 1 ); %#ok
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
    if dualized || ~dualize || pstr.dualize < 0,
        break
    end
    rsv(2:end) = 0;
    [ ncur, mcur ] = size(dbcA);
    dbcAT = dbcA';
    [ rows, cols ] = cvx_eliminate_mex( dbcAT, 1, ineqs, rsv ); %#ok
    if pstr.dualize > 0 || mcur - nrsv > ncur - nnz(rows),
        dbcA = -dbcAT;
        dir  = -dir;
        tmp  = Q; Q = P; P = tmp;
        tmp  = ineqs; ineqs = rsv; rsv = tmp;
        nnew = mcur;
        dualized = true;
        n = nnew - nrsv;
        for k = 1 : length(cones),
            [nn,nv] = size(cones(k).indices);
            nnnv = nn * nv;
            ndxs = n + 1 : n + nnnv;
            cones(k).indices = reshape(ndxs,nn,nv); %#ok
            n = n + nnnv;
        end
        ndxs = 1 : nnew;
    else
        clear dbcAT
        dbcA = dbcA(:,1:mcur-nrsv);
        P = P(:,1:mcur-nrsv);
        break;
    end
end

%
% STEP 4: Look for redundant slack variables. Look for situations like
%   minimize ( ... + x )
%   s.t. x >= 0
% or situations like
%   a' * x + s1 + s2 <= b
%   s1, s2 >= 0
% And eliminate the redundancy analytically.
%

nfin = length(ndxs);
ndxi = zeros(nnew,1);
ndxi(ndxs) = 1 : nfin;
if u_nneg,
    slacks = false(nfin,1);
    nneg = false(nfin,1);
end
do_slack = false;
for k = 1 : ncones,
    temp = cones(k).indices;
    temp = reshape(ndxi(temp),size(temp));
    cones(k).indices = temp; %#ok
    if u_nneg,
        temp = temp(cones(k).slacks>0,:);
        tmpv = sum(dbcA(temp,:)~=0,2)==1;
        if isequal(cones(k).type,'nonnegative'),
            tmpv = tmpv & (dbcA(temp,1)<=0);
            if any(tmpv), 
                do_slack = true; 
                slacks(temp) = tmpv;
                nneg(temp) = tmpv; 
            end
        else
            tmpv = tmpv & (dbcA(temp,1)==0);
            if any(tmpv),
                slacks(temp) = tmpv;
            end
        end
    end
end
if do_slack,
    t_slack = dbcA( slacks, : );
    d_slack = sum(t_slack<0,1);
    c_slack = max(1,sum(t_slack>0,1))-max(1,d_slack);
    c_slack(1) = -d_slack(1);
    if any( c_slack ),
        cndxs = c_slack ~= 0;
        t_slack = dbcA(nneg,cndxs);
        if nnz( t_slack ),
            [ rx, cx ] = find( max( 0, bsxfun( @times, t_slack, c_slack(:,cndxs) ) ) );
            if ~isempty( rx ),
                rndxs = find(nneg);
                cndxs = find(cndxs);
                cx = cndxs(cx); cx = cx(:);
                rx = rndxs(rx);
                dx = diff([0;find(diff(cx))]);
                dx = cumsum(1+sparse(1,cumsum(dx)+1,-dx,1,length(cx)));
                rx = rx(dx<=abs(c_slack(cx)));
                Q(:,rx) = [];
                dbcA(rx,:) = [];
                ndxi = ones(nfin,1);
                ndxi(rx) = 0;
                ndxi = cumsum(ndxi);
                ndxi(rx) = 0;
                cone2 = [];
                for k = 1 : ncones,
                    cone  = cones(k);
                    ctype = cone.type;
                    ndxs  = cone.indices;
                    ndxs  = reshape(ndxi(ndxs),size(ndxs));
                    if isequal(ctype,'nonnegative'),
                        ndxs = nonzeros(ndxs)';
                        if isempty(ndxs), continue; end
                    end
                    cone2 = cvx_pushcone( cone2, ctype, ndxs );
                end
                cones = cone2;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERSION / COMPATIBILITY CHECK %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~nargin || isempty( config )
    return
end
Pr = []; Pc = []; Pv = [];
rejected = {};
if nargin && ~isempty( config ),
    for k = 1 : length(cones),
        cone = cones(k);
        ctyp = cone.type;
        if isfield( config, ctyp ), 
            cmode = config.(ctyp);
            if cmode > 0, continue; end
        else
            cmode = 0;
        end
        ndxs = cone.indices;
        [nn,nv] = size(ndxs);
        success = false;
        switch ctyp,
            case 'absolute_value',
                if isfield( config, 'nonnegative' ) && config.nonnegative,
                    beta = sqrt(0.5);
                    ndx1 = ndxs(1,:);
                    ndx2 = ndxs(2,:);
                    Pr = [Pr,ndx1,ndx1,ndx2,ndx2]; %#ok
                    Pc = [Pc,ndx1,ndx2,ndx1,ndx2]; %#ok
                    Pv = [Pv,beta*[-ones(1,nv),ones(1,3*nv)]]; %#ok
                    ctyp = 'nonnegative';
                    ndxs = ndxs(:)';
                    success = true;
                end
            case 'rotated_lorentz',
                if isfield( config, 'lorentz' ) && config.lorentz,
                    % Convert to lorentz
                    ndx1 = ndxs(end-1,:);
                    ndx2 = ndxs(end,:);
                    beta = sqrt(0.5);
                    Pr = [Pr,ndx1,ndx1,ndx2,ndx2]; %#ok
                    Pc = [Pc,ndx1,ndx2,ndx1,ndx2]; %#ok
                    Pv = [Pv,beta*[-ones(1,nv),ones(1,3*nv)]]; %#ok
                    ctyp = 'lorentz';
                    success = true;
                end
            case 'semidefinite',
                if size( ndxs, 1 ) <= 3,
                    % Convert to rotated_lorentz or lorentz
                    ndx1 = ndxs(1,:);
                    ndx2 = ndxs(2,:);
                    ndx3 = ndxs(3,:);
                    if isfield( config, 'rotated_lorentz' ) && config.rotated_lorentz,
                        beta = sqrt(0.5);
                        Pr = [Pr,ndx1,ndx2,ndx3]; %#ok
                        Pc = [Pc,ndx2,ndx1,ndx3]; %#ok
                        Pv = [Pv,ones(1,nv)/beta,ones(1,nv),ones(1,nv)/beta]; %#ok
                        ctyp = 'rotated_lorentz';
                        success = true;
                    elseif isfield( config, 'lorentz' ) && config.lorentz,
                        Pr = [Pr,ndx1,ndx1,ndx2,ndx3,ndx3]; %#ok
                        Pc = [Pc,ndx2,ndx3,ndx1,ndx2,ndx3]; %#ok
                        Pv = [Pv,-ones(1,nv),ones(1,4*nv)]; %#ok
                        ctyp = 'lorentz';
                        success = true;
                    end
                end
            case 'hermitian_semidefinite',
                if size( ndxs, 1 ) <= 4,
                    ndx1 = ndxs(1,:);
                    ndx2 = ndxs(2,:);
                    ndx3 = ndxs(3,:);
                    ndx4 = ndxs(4,:);
                    if isfield( config, 'rotated_lorentz' ) && config.rotated_lorentz,
                        beta = sqrt(0.5);
                        Pr = [Pr,ndx1,ndx2,ndx3,ndx4]; %#ok
                        Pc = [Pc,ndx3,ndx1,ndx2,ndx4]; %#ok
                        Pv = [Pv,ones(1,nv)/beta,ones(1,2*nv),ones(1,nv)/beta]; %#ok
                        ctyp = 'rotated_lorentz';
                        success = true;
                    elseif isfield( config, 'lorentz' ) && config.lorentz,
                        Pr = [Pr,ndx1,ndx1,ndx2,ndx3,ndx4,ndx4]; %#ok
                        Pc = [Pc,ndx3,ndx4,ndx1,ndx2,ndx3,ndx4]; %#ok
                        Pv = [Pv,-ones(1,nv),ones(1,5*nv)]; %#ok
                        ctyp = 'lorentz';
                        success = true;
                    end
                end
                if ~success && isfield( config, 'semidefinite' ) && config.semidefinite, 
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
                    ctyp = 'semidefinite';
                    nfin = ndxs(end);
                    success = true;
                end
            case 'exponential',
                if dualize && isfield( config, 'lorentz' ),
                    success = true;
                end
        end
        if success,
            cones(k).type = ctyp; %#ok
            cones(k).indices = ndxs; %#ok
        elseif ~cmode,
            temp = cone.type;
            temp(temp=='_') = ' ';
            ordr = size( ndxs, 1 );
            if ordr == 1, 
                rejected{end+1} = sprintf( '%s variables', temp ); %#ok
            else
                rejected{end+1} = sprintf( 'order-%d %s cones', ordr, temp ); %#ok
            end
        end
    end
end
if ~isempty( rejected ),
    if nargin < 2 || isempty( name ),
       name = 'The current solver';
    end
    if length( rejected ) == 1,
        cvx_throw( 'CVX:IncompatibleSolver:%s does not support %s.\nPlease select a different solver.', name, rejected{1} );
    else
        temp = sprintf('\n    %s', rejected{:} );
        cvx_throw( 'CVX:IncompatibleSolver:%s does not support the following nonlinearities:%s\nPlease select a different solver.', name, temp );
    end
elseif u_exp && u_int,
    cvx_throw( 'Exponential variables and integer variables cannot both be used.' );
end
if ~isempty(Pr),
    Pd = true(1,size(dbcA,1),1);
    Pd(Pr) = false;
    Pd = find(Pd);
    Pr = [Pr,Pd];
    Pc = [Pc,Pd];
    Pv = [Pv,ones(size(Pd))];
    Pd = sparse(Pr,Pc,Pv);
    dbcA = Pd' * dbcA;
    Q = Q * Pd;
end
            
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
