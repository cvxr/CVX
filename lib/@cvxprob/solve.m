function solve( prob )

global cvx___
p = index( prob );
pr = cvx___.problems(p);
nobj = numel(pr);
if nobj > 1 & ~pr.separable,
    error( 'Non-separable multiobjective problems are not supported.' );
end
quiet = cvx___.problems(p).quiet;
obj = cvx___.problems(p).objective;
gobj = cvx___.problems(p).geometric;
clear pr
[ At, cones, sgn, Q, P, dualized ] = eliminate( prob, true );
% if any( strcmp( { cones.type }, 'exponential' ) ),
%     error( 'Constraints involving exp() and log() are not yet supported.' );
% end

dbca = At;
c = At( :, 1 );
At( :, 1 ) = [];
d = c( 1, : );
c( 1, : ) = [];
[ n1, m ] = size( At );
if n1 < 1,
    b = zeros( m, 1 );
else
    b = - At( 1, : ).';
    At( 1, : ) = [];
end
n = n1 - 1;
for k = 1 : length( cones ),
    cones(k).indices = cones(k).indices - 1;
end

%
% Ferret out the degenerate and overdetermined problems
%

tt = ( b' ~= 0 ) & ~any( At, 1 );
infeas = any( tt );
if m > n & n > 0,
    
    %
    % Overdetermined problem
    %
    
    x      = NaN * ones( n, 1 );
    y      = NaN * ones( m, 1 );
    oval   = NaN;
    status = 'Overdetermined';
    estr = sprintf( 'Overdetermined equality constraints detected.\n   CVX cannot solve this problem; but it is likely infeasible.' );
    if ~quiet,
        disp( estr );
    else
        warning( estr );
    end
    pval = NaN;
    dval = NaN;

elseif n ~= 0 & ~infeas & ( any( b ) | any( c ) ),
        
    %
    % Call solver
    %
    
    prob = cvx___.problems( p );
    solv = prob.solver;
    lsolv = lower(solv);
    prec = prob.precision;
    sfunc  = [ 'cvx_solve_', lsolv ];
    tt = find( strcmp( { cones.type }, 'exponential' ) );
    need_iter = ~isempty( tt );
    if ~quiet,
        disp( ' ' );
        spacer = '-';
        if need_iter,
            disp( 'Successive approximation method to be employed.' );
            disp( sprintf( '   %s will be called several times to refine the solution.', solv ) );
            disp( sprintf( '   Original size: %d variables, %d equality constraints', n, m ) );
            spacer = spacer(:,ones(1,65));
        else
            disp( sprintf( 'Calling %s: %d variables, %d equality constraints', solv, n, m ) );
            spacer = spacer(:,ones(1,60));
        end
        if dualized,
            disp( sprintf( '   For improved efficiency, %s is solving the dual problem.', solv ) );
        end
        if ~need_iter,
            disp( spacer );
        end
    end
    cvx_setspath( solv );
    if cvx___.profile, profile off; end
    if need_iter,
        
        %
        % Cone:
        %     cl { (x,y,z) | y*exp(x/y) <= z, y > 0 }
        %   = cl { (x,y,z) | x <= -y*log(y/z), z > 0 }
        % Approximation: given a shift point x0,
        %    { (x,y,z) | y*exp(x0)*pos(1+(x/y-x0)/16)^16 <= z, y > 0 }
        %    { (x,y,z) | y+(x-x0*y)/16 <= exp(-x0/16)*geomean([z,y],[],[1,15])
        % Transformed cone:
        %   4 semidefinite cones, 1 lorentz cone, 1 slack
        %   [ w1    ][ w4    ] [ w7    ] [ w10     ] [ w13 ]
        %   [ w2 w3 ][ w5 w6 ] [ w8 w9 ] [ w11 w12 ] [ w14 ] w15
        %   w2 = w4, w5 = w7, w8 = w10
        %   w3 = w6, w6 = w9, w9 = w12, w12 = w14 / xw
        %   exp(-x0/16) * w11 = w3 + w13 / 16 + w15
        % Recovery:
        %   x = w13 + x0 * x3
        %   y = w3
        %   z = w1
        %
        
        ndxs  = cat( 2, cones(tt).indices );
        cones(tt) = [];
        nc    = size(ndxs,2);
        xndxs = ndxs(1,:);
        yndxs = ndxs(2,:);
        zndxs = ndxs(3,:);
        x0    = realmin * ones(nc,1);
        maxw  = log(realmax);
        
        epow = 8;
        switch epow,
            case 16,
                QAi  = [ 2, 4, 3, 6, 5, 7, 6, 9, 8,10, 9,12,12,    14, 3,    11,  13, 15 ]';
                QAj  = [ 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7,     7, 8,     8,   8,  8 ]';
                QAv  = [+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-0.123, 1,-0.234,1/16,  1 ]';
                QAr  = [3,4,2,5,6,7,8,9,10,11,12,13,1,14,15];
                ewid = 1.75;
            case 8,
                QAi  = [ 2, 4, 3, 6, 5, 7, 6, 9, 9,    11, 3,     8,  10, 12 ]';
                QAj  = [ 1, 1, 2, 2, 3, 3, 4, 4, 5,     5, 6,     6,   6,  6 ]';
                QAv  = [+1,-1,+1,-1,+1,-1,+1,-1,+1,-0.123, 1,-0.234,1/8,   1 ]';
                QAr  = [3,4,2,5,6,7,8,9,10,1,11,12];
                ewid = 1.22;
            case 4,
                QAi  = [ 2, 4, 3, 6, 6,     8, 3,     5,   7,  9 ]';
                QAj  = [ 1, 1, 2, 2, 3,     3, 4,     4,   4,  4 ]';
                QAv  = [+1,-1,+1,-1,+1,-0.123, 1,-0.234, 1/4,  1 ]';
                QAr  = [3,4,2,5,6,7,1,8,9];
                ewid = 0.84;
        end
        
        nQA     = max(QAi);
        mQA     = max(QAj);
        nc      = size(ndxs,2);
        new_n   = n + (nQA-3) * nc;
        new_m   = m + mQA * nc;
        n_ndxs  = [ ndxs ; reshape( n + 1 : new_n, nQA-3, nc  ) ];
        n_ndxs  = n_ndxs(QAr,:);
        if ~quiet,
            disp( sprintf( '   Approximation size: %d variables, %d equality constraints', new_n, new_m ) );
            disp( spacer );
        end
        
        % Perform (x,y,z) ==> w transformation on A and C
        c (new_n,end) = 0;
        At(new_n,end) = 0;
        b (new_m,end) = 0;
        
        % Add new cone constraints
        lQA = length(QAi);
        nc0 = 0:mQA:mQA*(nc-1);
        nc1 = ones(1,nc);
        At = [ At, sparse( n_ndxs(QAi,:), ...
            QAj(:,nc1) + nc0(ones(lQA,1),:), ...
            QAv(:,nc1), new_n, mQA * nc ) ];
        
        wndxs = n_ndxs(nQA-1,:) + (m+mQA-2+nc0) * new_n;
        endxs = n_ndxs(nQA-4,:) + (m+mQA-1+nc0) * new_n;
        yorig_A = At(yndxs,:);
        yorig_C = c(yndxs,:);
        xorig_A = At(xndxs,:);
        xorig_C = c(xndxs,:);
        xorig_A(:,m+1:end) = 0;
        
        ncone.type       = 'semidefinite';
        ncone.indices    = reshape(n_ndxs(1:nQA-3,:),3,(nQA-3)*nc/3);
        ncone(2).type    = 'lorentz';
        ncone(2).indices = reshape(n_ndxs(nQA-2:nQA-1,:),2,nc);
        ncone(3).type    = 'nonnegative';
        ncone(3).indices = n_ndxs(nQA,:);
        cones = [ cones, ncone ];

        xw = epow;
        rel_err = 1e-6;
        last_err = Inf;
        
        prec_list = prec(3); % max( prec(3), 1e-2 );
        for k = 3: -1 : 1,
            if prec(k) > 0,
                if 0, % prec(k) < 0.1 * prec_list(end),
                    prec_list = [ prec_list, sqrt(prec(k)*prec_list(1)), prec(k) ];
                elseif prec(k) < prec_list(end),
                    prec_list = [ prec_list, prec(k) ];
                end
            end
        end
        if prec(1) == 0,
            % prec_list(end+1) = prec_list(end) / 8;
            prec_list(end+1) = 0;
        end
        prec_ndx = 1;
        eshift   = 0;
        ninfeas  = 0;
        best_x   = [];
        best_px  = Inf;
        best_ox  = Inf;
        best_py  = Inf;
        best_oy  = -Inf;
        best_y   = [];
        solv_warn = false;
        use_init = false;
        nprec = prec_list(1) * [1,1,1];
        if ~quiet,
            disp( ' Target     Forcing     Conic    Solver' );
            disp( 'Precision    Bias       Error    Status' );
            disp( '---------------------------------------' );
        end
        best_x = NaN * ones(n,1);
        best_y = NaN * ones(m,1);
        dobj = -Inf;
        pobj = +Inf;
        stagnant = false;
        last_slow = false;
        dscale = blkdiag(speye(m,m),speye(new_m-m,new_m-m));
        for iter = 1 : 1000,
            dx0 = diag(sparse(x0));
            At(yndxs,:) = yorig_A + dx0 * xorig_A;
            c (yndxs,:) = yorig_C + dx0 * xorig_C;
            At(endxs) = -exp(-(x0+eshift)/epow);
            At(wndxs) = -1./xw;
            [ x, y, status ] = feval( sfunc, At, b, c, cones, true, nprec );
            x_good = ~isnan( x(1) );
            y_good = ~isnan( y(1) );
            tndx   = 2 + ( status(3) == 'a' );
            tprec  = nprec(tndx);
            if x_good,
                xxx = x(xndxs,:);
                yyy = x(yndxs,:);
                zzz = x(zndxs,:);
                nx0 = xxx ./ yyy;
                nx1 = log( zzz ./ yyy ) - x0;
                nxe = nx0 - nx1;
                x_clean = all( nxe < 0 );
                pob = c' * x;
            elseif any( eshift > 0 ),
                eshift = min( eshift, 0 );
                x_clean = false;
                nxe = 0;
            else
                nxe = 0;
                x_clean = true;
                pob = Inf;
            end
            if y_good,
                y(m+1:end) = 0;
                z = c - At * y;
                uuu = z( xndxs, : );
                vvv = z( yndxs, : );
                www = z( zndxs, : );
                nx2 = 1 - vvv ./ uuu;
                nx3 = log( - uuu ./ www ) - x0;
                nxe = nx0 - nx1;
                nxd = nx2 - nx3;
                y_clean = all( nxd > 0 ) | ~any( eshift );
                dob = b' * y;
                if x_good, % Solved
                    dx0 = 0.5 * ( nx0 + nx2 );
                else       % Infeasible
                    dx0 = 0.5 * ( nx2 + nx3 );
                end
            elseif any( eshift < 0 ),
                eshift = max( eshift, 0 );
                y_clean = false;
                nxd = 0;
            else
                nxd = 0;
                y_clean = true;
                dob = -Inf;
                if x_good, % Unbounded
                    dx0 = 0.5 * ( nx0 + nx1 );
                else       % Failed
                   break;
                end
            end
            if x_clean & ( tprec < best_px | ( tprec == best_px & pob < best_ox ) ),
                best_x  = x(1:n);
                best_x(xndxs,:) = xxx + x0 .* yyy;
                best_px = tprec;
                best_ox = pob;
            end
            if y_clean & ( tprec < best_py | ( tprec == best_py & dob > best_oy ) ),
                best_y  = y(1:m);
                best_py = tprec;
                best_oy = dob;
            end
            err = 0;
            if best_px > tprec,
                err = max( err, max(nxe) );
            end
            if best_py > tprec,
                err = max( err, -min(nxd) );
            end
            slow = err > 0.9 * last_err;
            stagnant = stagnant | ( slow & last_slow );
            if ~quiet,
                if stagnant, stagc = 'S'; else stagc = ' '; end
                disp( sprintf( '%9.3e  %9.3e%c %9.3e  %s', nprec(1), max(eshift), stagc, err, status ) );
            end
            if ( best_px == tprec & best_py == tprec ) | nprec(1) > prec(3),
                if nprec(1) == prec(1) | tndx == 3,
                    break;
                end
                prec_ndx = prec_ndx + 1;
                nprec(1) = prec_list(prec_ndx);
                nprec(2) = max( nprec(1), prec(2) );
                nprec(3) = max( nprec(1), prec(3) );
                stagnant = false;
                err = Inf;
                eshift = 0;
            elseif stagnant & best_px > tprec,
                eshift = eshift + 2 * max(nxe,0);
            elseif stagnant & best_py > tprec,
                eshift = 0;
            elseif slow,
                eshift = eshift + (nxe>0).*min(nxd,nxe) + (nxd<0).*max(nxd,nxe);
            end
            last_err = err;
            last_slow = slow;
            xw = max( abs( dx0 ), ewid );
            x0 = x0 + max( dx0, -epow );
        end
        if isnan( best_x(1) ), 
            status = 'Infeasible',
        elseif isnan( best_y(1) ), 
            status = 'Unbounded';
        else
            status = 'Solved';
        end
        best_p = max( best_px, best_py );
        if best_p > prec(3),
            status = 'Failed';
        elseif best_p > prec(2),
            status = [ 'Inaccurate/', status ];
        end
        x = best_x;
        y = best_y;
        c = c(1:n,:);
        c(yndxs,:) = yorig_C;
    else
        [ x, y, status ] = feval( sfunc, At, b, c, cones, quiet, prec );
    end
    if cvx___.profile, profile on; end
    if ~cvx___.path.hold, 
        cvx_setspath(''); 
    end
    switch status,
    case { 'Solved', 'Inaccurate/Solved' },
        oval = sgn * ( c' * x + d' );
        pval = 1;
        dval = 1;
    case { 'Infeasible', 'Inaccurate/Infeasible' },
        oval = sgn * Inf;
        pval = NaN;
        dval = 0;
    case { 'Unbounded', 'Inaccurate/Unbounded' },
        oval = -sgn * Inf;
        pval = 0;
        dval = NaN;
    otherwise,
        oval = NaN;
        pval = NaN;
        dval = NaN;
    end
    if ~quiet,
        disp( spacer );
    end
    
elseif infeas,
    
    %
    % Infeasible
    %
    
    if ~quiet,
        disp( 'Trivial infeasibilities detected; solution determined analytically.' );
    end
    status = 'Infeasible';
    x = NaN * ones( n, 1 );
    b( ~tt ) = 0;
    y = - b / ( b' * b );
    oval = sgn * Inf;
    pval = NaN;
    dval = 0;
    
else
    
    %
    % The origin is optional
    %
    
    if ~quiet,
        disp( 'Homogeneous problem detected; solution determined analytically.' );
    end
    status = 'Solved';
    x = zeros( n, 1 );
    y = zeros( m, 1 );
    oval = sgn * d;
    pval = 1;
    dval = 1;
    
end

if dualized,
    switch status,
        case 'Infeasible', status = 'Unbounded';
        case 'Unbounded',  status = 'Infeasible';
        case 'Inaccurate/Infeasible', status = 'Inaccurate/Unbounded';
        case 'Inaccurate/Unbounded',  status = 'Inaccurate/Infeasible';
    end
end

if ~quiet,
    disp( sprintf( 'Status: %s', status ) );
end

cvx___.problems( p ).status = status;

%
% Push the results into the master CVX workspace
%

x = full( Q * [ pval ; x ] );
y = full( P * [ dval ; y ] );
if dualized,
    cvx___.x = y;
    cvx___.y = x(2:end);
else
    cvx___.x = x;
    cvx___.y = y(2:end);
end
if nnz( cvx___.exponential ),
    esrc = find( cvx___.exponential );
    edst = cvx___.exponential( esrc );
    cvx___.x( edst ) = exp( cvx___.x( esrc ) );
end

%
% Compute the objective
%

if ~isempty( obj ),
    if isinf( oval ) | isnan( oval ),
        oval = oval * ones(size(obj));
    else
        oval = cvx_value( obj );
    end
    oval(gobj) = exp(oval(gobj));
end
oval = full(oval);
cvx___.problems( p ).result = oval;
if ~quiet,
    if length( oval ) == 1,
        disp( sprintf( 'Optimal value (cvx_optval): %+g', oval ) );
    else
        disp( sprintf( 'Optimal value (cvx_optval): (multiobjective)' ) );
    end
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
