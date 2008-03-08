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
    spacer = '-';
    spacer = spacer(:,ones(1,60));
    sfunc  = [ 'cvx_solve_', lsolv ];
    tt = find( strcmp( { cones.type }, 'exponential' ) );
    need_iter = ~isempty( tt );
    if ~quiet,
        disp( ' ' );
        disp( sprintf( 'Calling %s: %d variables, %d equality constraints', solv, n, m ) );
        if dualized,
            disp( sprintf( 'Note: for improved efficiency, %s is solving the dual problem.', solv ) );
        end
        if need_iter,
            disp( sprintf( 'Note: Successive approximation method to be employed;\n   %s will be called several times to refine the solution.', solv ) );
        end
        disp( spacer );
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
            disp( sprintf( 'Approximation size: %d variables, %d equality constraints', new_n, new_m ) );
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
        
        prec_list = max( prec(3), 1e-2 );
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
            prec_list(end+1) = prec_list(end) / 8;
            prec_list(end+1) = 0;
        end
        prec_ndx = 1;
        best_ndx = 0;
        eshift   = 0;
        ninfeas  = 0;
        best_x   = [];
        best_y   = [];
        best_s   = [];
        failed_flag = false;
        prec_bump = false;
        solv_warn = false;
        use_init = false;
        nprec = prec_list(1) * [1,1,1];
        if ~quiet,
            disp(sprintf('Target precision: %e', nprec(1)));
        end
        XYZ = {};
        for iter = 1 : 1000,
            if prec_bump,
                prec_bump = false;
            else
                dx0 = diag(sparse(x0));
                At(yndxs,:) = yorig_A + dx0 * xorig_A;
                c (yndxs,:) = yorig_C + dx0 * xorig_C;
                At(endxs) = -exp(-(x0+eshift)/epow);
                At(wndxs) = -1./xw;
            end
            XYZ = {};
            use_init = ~isempty( XYZ );
            [ x, y, status ] = feval( sfunc, At, b, c, cones, quiet, nprec );
            % [ x, y, status, XYZ ] = feval( sfunc, At, b, c, cones, quiet, nprec, XYZ );
            xxx = x(xndxs,:);
            yyy = x(yndxs,:);
            switch status,
                case { 'Solved', 'Inaccurate/Solved', 'Unbounded', 'Inaccurate/Unbounded' },
                    ninfeas = 0;
                    nx0 = xxx ./ yyy;
                    nx1 = log( x(zndxs,:) ./ yyy ) - x0;
                    nxe = max( 0, nx0 - nx1 );
                    dx0 = 0.5 * ( nx0 + nx1 );
                    if nprec(1) == 0,
                        sms = nxe ~= 0;
                    else
                        sms = abs( dx0 ) < 0.05 * nxe | eshift;
                    end
                    err = max( nxe );
                    if ~quiet,
                        disp(sprintf('Approxmation error: %e',err));
                    end
                    if err == 0 | nprec(1) > prec(2),
                        if err == 0,
                            best_y = y;
                            best_x = x;
                            best_x(xndxs,:) = xxx + x0 .* yyy;
                            best_s = status;
                            best_ndx = prec_ndx;
                            if nprec(1) == prec(1), 
                                break; 
                            end
                        end
                        prec_ndx = prec_ndx + 1;
                        nprec(1) = prec_list(prec_ndx);
                        nprec(2:3) = max( nprec(1), prec(2) );
                        last_err = Inf;
                        if ~quiet,
                            disp(sprintf('Target precision: %e',nprec(1)));
                        end
                        if any( eshift ),
                            eshift = 0;
                            XYZ = {};
                        elseif 0, % err == 0,
                            prec_bump = true;
                            continue;
                        end
                    elseif any( sms ),
                        if ~solv_warn,
                            warning( ...
                                sprintf( ...
[ 'Your problem has activated a weakly tested section of the successive\n',...
  'approximation code. Please let Michael Grant know that you have seen\n',...
  'this message, and if possible provide him with the model that produced it.' ] ) );
                            solv_warn = true;
                        end
                        eshift = eshift + nxe .* sms;
                        last_err = 0;
                    else
                        last_err = err;
                    end
                    dx0 = 0.5 * ( nx0 + nx1 ); % max( min( nx0, nx1 ), 0 ) + min( max( nx0, nx1 ), 0 );
                    xw  = max( ewid, 2 * abs( nx0 - nx1 ) );
                    x0  = x0 + max( dx0, -epow );
                case { 'Infeasible', 'Inaccurate/Infeasible' },
                    ninfeas = ninfas + 1;
                    if any( eshift ),
                        eshift = eshift * ( 0.5 * ninfeas < 3 );
                    elseif norm([x0-xw;x0+xw],Inf) > maxw,
                        best_y = y;
                        best_x = x;
                        best_s = status;
                        best_ndx = length(prec_list);
                        break;
                    end
                    if ~quiet,
                        disp(sprintf('Infeasible approximation; widening trust region'));
                    end
                    xw = xw * 8;
                case 'Failed',
                    if use_init,
                        XYZ = {};
                    else
                        if best_ndx == 0,
                            best_x = x;
                            best_y = y;
                        end
                        break;
                    end
            end
        end
        if best_ndx == 0 | prec_list(best_ndx) > prec(3),
            status = 'Failed';
        elseif prec_list(best_ndx) > prec(2),
            status = [ 'Inaccurate/', best_s ];
        elseif strcmp( best_s(1:3), 'Ina' ),
            status = best_s(12:end);
        else
            status = best_s;
        end
        x = best_x(1:n,:);
        y = best_y(1:m,:);
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

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
