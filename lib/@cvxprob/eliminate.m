function [ dbCA, cones, dir, Q, P, exps, dualized ] = eliminate( prob, destructive, config )
if nargin < 2, destructive = false; end

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

[ dbCA, cones, dir, Q, P, exps ] = extract( prob, destructive );
dualized = false;
if size( dbCA, 1 ) == 1, 
    return; 
end

%
% Negate the objective so that the transformation matrices P and Q are
% properly formed.
%

can_dual = config.dualize;
dbCA(:,1) = -dbCA(:,1);
for pass = 1  : 2,
    
    nold = size(dbCA,1);
    rsv  = zeros( nold, 1 );
    nneg = false( nold, 1 );
    slacks = false( nold, 1 );
    for k = 1 : length( cones ),
        ctyp = cones(k).type;
        ndxs = cones(k).indices;
        rsv(ndxs) = 1;
        switch ctyp,
            case 'nonnegative',
                slacks(ndxs) = true;
                nneg(ndxs) = true;
            case { 'exponential', 'lorentz' },
                slacks(ndxs(end,:)) = true;
            case 'rotated-lorentz',
                slacks(ndxs(end-1:end,:)) = true;
            case 'semidefinite',
                q = round(0.5*(sqrt(8*size(ndxs,1)+1)-1));
                slacks(ndxs(cumsum([1,q:-1:2]),:)) = true;
            case 'hermitian-semidefinite',
                q = round(sqrt(size(ndxs,1)));
                slacks(ndxs(cumsum([1,2*q-1:-2:2]),:)) = true;
        end
    end
    ndxs = ( 1 : nold )';
    
    ineqs = zeros(1,size(dbCA,2));
    ineqs(1) = 1;
    rsv(1) = 1;

    last_success = 3;
    while true,
        
        %
        % STEP 1: Look for columns which differ only by a constant factor.
        % These correspond to redundant equality constraints. These occur
        % often enough as as consequence of our tranformation method, and
        % they cause problems in solvers, so we must eliminate them. Of
        % course, if there are more complex linear dependencies in the
        % equality constraints, we can't do anything about that.
        %
        
        [ xR, dbCA ] = cvx_bcompress( dbCA, 'full', 1 );
        if size( xR, 1 ) ~= size( xR, 2 ),
            last_success = 1;
            P       = P * cvx_invert_structure( xR );
            ineqs   = ( xR * ineqs(:) )' ~= 0;
            ineqs   = +ineqs;
        elseif last_success == 1,
            break;
        end
        
        %
        % STEP 2: Look for variables that we can eliminate without
        % increasing fill-in. This means looking for rows or columns
        % with only 1, 2, or (in some cases) 3 nonzeros.
        %

        success = false;
        while true,
            [ rows, cols ] = cvx_eliminate_mex( dbCA, 1, rsv, ineqs );
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
            A11  = dbCA( rowX, colX );
            A12  = dbCA( rowX, cols );
            A21  = dbCA( rows, colX );
            A22  = dbCA( rows, cols );
            if ( size( A22, 1 ) ~= size( A22, 2 ) || nnz( A22 ) ~= size( A22, 1 ) ),
                error( 'There seems to be an error in the CVX presolver routine.\nPlease report this to the authors; and if possible, include the\ncvx model and data that gave you this error.', 1 ); %#ok
            end
            [ ii, jj, vv ] = find( A22 );
            A22i   = sparse( jj, ii, 1.0 ./ vv );
            temp   = - A22i * A21;
            P      = P( :, colX ) + P( :, cols ) * temp;
            temp   = - A12 * A22i;
            Q      = Q( :, rowX ) + Q( :, rows ) * temp';
            dbCA   = A11 + temp * A21;
            rsv    =   rsv( rowX, : );
            ndxs   =  ndxs( rowX, : );
            slacks = slacks( rowX, : );
            nneg   = nneg( rowX, : );
        end
        if success,
            last_success = 2;
        elseif last_success == 2,
            break;
        end
        
        %
        % STEP 3: Look for redundant slack variables.
        %
        
        success = false;
        t_slack = bsxfun( @times, dbCA, slacks );
        t_slack = bsxfun( @times, t_slack, ( sum( t_slack ~= 0, 2 ) == 1 ) & ( t_slack(:,1) == 0 ) );
        if nnz( t_slack ),
            c_slack = max(0,sum(t_slack>0,1)-1)-max(0,sum(t_slack<0,1)-1);
            c_slack(1) = 0; 
            if nnz( c_slack ),
                [ rx, cx ] = find( max( 0, bsxfun( @times, nneg, bsxfun( @times, t_slack, c_slack ) ) ) );
                if ~isempty( rx ),
                    dx = diff([0;find(diff(cx))]);
                    dx = cumsum(1+sparse(1,cumsum(dx)+1,-dx,1,length(cx)));
                    rx = rx(dx<=abs(c_slack(cx)));
                    Q(:,rx) = [];
                    dbCA(rx,:) = [];
                    rsv(rx) = [];
                    ndxs(rx) = [];
                    nneg(rx) = [];
                    slacks(rx) = [];
                    success = true;
                end
            end
        end
        if success,
            last_success = 3;
        elseif last_success == 3,
            break;
        end
        
    end
    
    if ~can_dual || isempty(cones),
        break;
    end
    
    %
    % Check to see if dualization will result in smaller problem
    %
    
    rsv(1) = 0;
    n_ineq = nnz(any(dbCA(rsv&(sum(dbCA~=0,2)==1)&~dbCA(:,1),:),1));
    rsv(1) = 1;
    [n1,m1] = size(dbCA);
    m_pri = m1 - 1;
    n_pri = n1 - 1;
    m_dua = n1 - n_ineq - 1;
    n_dua = nnz( rsv ) + m1 - n_ineq - 1;
    if ( ( m_pri > n_pri ) || ( m_pri * m_pri * n_pri > m_dua * m_dua * n_dua ) ) && ( m_dua <= n_dua ),
        ndxi = zeros(nold,1);
        ndxi(ndxs) = 1 : length(ndxs);
        PP = cell(2,length(cones));
        n_cur = m1;
        prune = false;
        tt = true(1,length(cones));
        for k = 1 : length(cones),
            temp = cones(k).indices;
            temp = reshape(ndxi(temp),size(temp));
            if ~all(temp),
                if ~any(temp),
                    tt(k) = false;
                    prune = true;
                    continue;
                end
                temp = nonzeros(temp)';
            end
            [nn,nv] = size(temp);
            switch cones(k).type,
                case 'semidefinite',
                    nt = 0.5*(sqrt(8*nn+1)-1);
                    SS = 'symmetric';
                case 'hermitian-semidefinite',
                    nt = sqrt(nn);
                    SS = 'hermitian';
                case 'exponential',
                    persistent SS_exp %#ok
                    if isempty( SS_exp ),
                        SS_exp = sparse([0,-1,0;-1,0,0;0,0,exp(-1)]);
                    end
                    SS = cvx_replicate_structure(SS_exp,nv);
                otherwise,
                    SS = [];
            end
            PP{k} = sparse(1:numel(temp),max(temp,1),temp~=0,numel(temp),n1);
            if ~isempty(SS),
                if ischar(SS),
                    SS = cvx_create_structure([nt,nt,nv],SS);
                    SS = SS * SS';
                end
                PP{k} = SS * PP{k};
            end
            cones(k).indices = reshape(n_cur+1:n_cur+nn*nv,nn,nv);
            n_cur = cones(k).indices(end);
        end
        if prune,
            cones = cones(tt);
            PP = PP(tt);
        end
        dbCA  = vertcat(dbCA',PP{:});
        dir   = -dir;
        tmp   = Q; Q = P; P = tmp;
        nold  = size(dbCA,1);
        Q(:,nold) = 0;
        dualized = true;
    else
        break;
    end
    
end

%
% Move the cone indices to their new locations and transform
%

nnew = size(dbCA,1);
ncone = length(cones);
ndxi = zeros(nold,1);
ndxi(ndxs) = 1 : length(ndxs);
prune = false;
for k = 1 : ncone,
    cone    = cones(k);
    ndxs    = cone.indices;
    [nn,nv] = size(ndxs);
    ndxs    = reshape( ndxi(ndxs), nn, nv );
    if ~all( ndxs(:) )
        if k > 1 || ~isequal( cone.type, 'nonnegative' ),
            error( 'CVX:InternalError', 'Internal CVX presolve error. Please report this to CVX Support.' );
        end
        ndxs = nonzeros(ndxs)';
        prune = isempty( ndxs );
    end
    cone.indices = ndxs;
    cones(k) = cone;
end
if prune,
    cones = cones(2:end);
end
Pr = []; Pc = []; Pv = [];
for k = 1 : ncone,
    cone    = cones(k);
    ndxs    = cone.indices;
    switch cone.type,
        case 'nonnegative',
        case 'lorentz',
            if nn == 2 && config.preferLP,
                beta = sqrt(0.5);
                Pr = [Pr,ndxs(1,:),ndxs(1,:),ndxs(2,:),ndxs(2,:)]; %#ok
                Pc = [Pc,ndxs(1,:),ndxs(2,:),ndxs(1,:),ndxs(2,:)]; %#ok
                Pv = [Pv,beta*[-ones(1,nv),ones(1,3*nv)]]; %#ok
                cone.indices = cone.indices(:).';
                cone.type = 'nonnegative';
            end
        case 'rotated-lorentz',
            if config.unrotateSOCP,
                ndxs = ndxs(end-1:end,:);
                beta = sqrt(0.5);
                Pr = [Pr,ndxs(1,:),ndxs(1,:),ndxs(2,:),ndxs(2,:)]; %#ok
                Pc = [Pc,ndxs(1,:),ndxs(2,:),ndxs(1,:),ndxs(2,:)]; %#ok
                Pv = [Pv,beta*[-ones(1,nv),ones(1,3*nv)]]; %#ok
                cone.type = 'lorentz';
            end
        case 'semidefinite',
            if nn == 3 && config.preferSOCP,
                if config.unrotateSOCP,
                    Pr = [Pr,ndxs(1,:),ndxs(1,:),ndxs(2,:),ndxs(3,:),ndxs(3,:)]; %#ok
                    Pc = [Pc,ndxs(2,:),ndxs(3,:),ndxs(1,:),ndxs(2,:),ndxs(3,:)]; %#ok
                    Pv = [Pv,-ones(1,nv),ones(1,4*nv)]; %#ok
                    cone.type = 'lorentz';
                else
                    beta = sqrt(0.5);
                    Pr = [Pr,ndxs(1,:),ndxs(2,:),ndxs(3,:)]; %#ok
                    Pc = [Pc,ndxs(2,:),ndxs(1,:),ndxs(3,:)]; %#ok
                    Pv = [Pv,ones(1,nv)/beta,ones(1,nv),ones(1,nv)/beta]; %#ok
                    cone.type = 'rotated-lorentz';
                end
            elseif ~config.capableSDP,
                error( 'CVX:IncompatibleSolver', 'This solver does not support semidefinite constraints.' );
            end
        case 'hermitian-semidefinite',
            if nn == 4 && config.preferSOCP,
                if config.unrotateSOCP,
                    Pr = [Pr,ndxs(1,:),ndxs(1,:),ndxs(2,:),ndxs(3,:),ndxs(4,:),ndxs(4,:)]; %#ok
                    Pc = [Pc,ndxs(3,:),ndxs(4,:),ndxs(1,:),ndxs(2,:),ndxs(3,:)]; %#ok
                    Pv = [Pv,-ones(1,nv),ones(1,4*nv)]; %#ok
                    cone.type = 'lorentz';
                else
                    beta = sqrt(0.5);
                    Pr = [Pr,ndxs(1,:),ndxs(2,:),ndxs(3,:),ndxs(4,:)]; %#ok
                    Pc = [Pc,ndxs(3,:),ndxs(1,:),ndxs(2,:),ndxs(4,:)]; %#ok
                    Pv = [Pv,ones(1,nv)/beta,ones(1,2*nv),ones(1,nv)/beta]; %#ok
                    cone.type = 'rotated-lorentz';
                end
            elseif ~config.capableSDP,
                error( 'CVX:IncompatibleSolver', 'This solver does not support semidefinite constraints.' );
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
                ndxs = [ndxs;reshape(nnew+1:nnew+n2*(n2+1)*nv,[],nv)]; %#ok
                Pr = [ Pr, vec( [ ndxs(rs,:)  ; ndxs(rs,:)  ; ndxs(is,:)  ; ndxs(is,:)  ] )' ]; %#ok
                Pc = [ Pc, vec( [ ndxs(rd1,:) ; ndxs(rd2,:) ; ndxs(id1,:) ; ndxs(id2,:) ] )' ]; %#ok
                Pv = [ Pv, vec( [ ones(2*numel(rs)+numel(is),nv) ; -ones(numel(is),nv) ] )' ]; %#ok
                cone.indices = ndxs;
                cone.type = 'semidefinite';
            end
        otherwise,
            if isequal( cone.type(1:2), 'i_' ) && ~config.INTcapable,
                error( 'CVX:IncompatibleSolver', 'This solver does not support integer constraints.' );
            end
    end
    cones(k) = cone;
end
if ~isempty(Pr),
    Pd = true(1,nnew);
    Pd(Pr) = false;
    Pd = find(Pd);
    Pr = [Pr,Pd];
    Pc = [Pc,Pd];
    Pv = [Pv,ones(size(Pd))];
    Pd = sparse(Pr,Pc,Pv);
    dbCA = Pd' * dbCA;
    Q    = Q * Pd;
    tt = strcmp( { cones.type }, 'nonnegative' );
    nt = nnz( tt );
    if nt > 1 || ( nt == 1 && ~isequal( cones(1).type, 'nonnegative' ) ),
        indices = cat(2,cones(tt).indices);
        cones = [ struct( 'type', 'nonnegative', 'indices', indices ), cones(~tt) ];
    end
end

%
% Return the objective back to normal.
%

if dualized,
    P = -P;
    P(:,1) = -P(:,1);
else
    dbCA(:,1) = -dbCA(:,1);
end
            
            
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
