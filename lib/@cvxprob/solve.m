function solve( prob )

global cvx___
p = index( prob );
pr = cvx___.problems(p);
nobj = numel(pr.objective);
if nobj > 1 & ~pr.separable,
    error( 'Non-separable multiobjective problems are not supported.' );
end
quiet = cvx___.problems(p).quiet;
obj   = cvx___.problems(p).objective;
gobj  = cvx___.problems(p).geometric;
clear pr
[ At, cones, sgn, Q, P, dualized ] = eliminate( prob, true );

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
    if dualized,
        status = 'Underdetermined';
        estr = sprintf( 'Underdetermined inequality constraints detected.\n   CVX cannot solve this problem; but it is likely unbounded.' );
    else
        status = 'Overdetermined';
        estr = sprintf( 'Redundant equality constraints detected.\n   CVX cannot solve this problem; but it is likely infeasible.' );
    end
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
    if isempty( cones ),
        texp = [];
    else
        texp = find( strcmp( { cones.type }, 'exponential' ) );
    end
    need_iter = ~isempty( texp );
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
        %    { (x,y,z) | y+(x-x0*y)/16 <= exp(-x0/16)*geo_mean([z,y],[],[1,15])
        % Transformed cone:
        %   4 semidefinite cones, 1 free, 1 slack
        %   [ w1    ][ w4    ] [ w7    ] [ w10     ] w13
        %   [ w2 w3 ][ w5 w6 ] [ w8 w9 ] [ w11 w12 ] w14
        %   w2 = w4, w5 = w7, w8 = w10
        %   w3 = w6, w6 = w9, w9 = w12,
        %   exp(-x0/16) * w11 = w3 ( 1 - x0 / 16 ) + w13 / 16 + w14
        % Recovery:
        %   x = w13
        %   y = w3
        %   z = w1
        %
        
        ndxs  = cat( 2, cones(texp).indices );
        nc    = size(ndxs,2);
        xndxs = ndxs(1,:);
        yndxs = ndxs(2,:);
        zndxs = ndxs(3,:);
        x0    = realmin * ones(nc,1);
        maxw  = log(realmax);
        
        epow = 8;
        switch epow,
            case 16,
                QAi  = [ 2, 4, 3, 6, 5, 7, 6, 9, 8,10, 9,12,3,        11,  13, 14 ]';
                QAj  = [ 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,7,         7,   7,  7 ]';
                QAv  = [+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,0.123,-0.234,1/16,  1 ]';
                QAr  = [3,4,2,5,6,7,8,9,10,11,12,13,1,14];
                ewid = 1.75;
            case 8,
                QAi  = [ 2, 4, 3, 6, 5, 7, 6, 9,3,         8,  10, 11 ]';
                QAj  = [ 1, 1, 2, 2, 3, 3, 4, 4,5,         5,   5,  5 ]';
                QAv  = [+1,-1,+1,-1,+1,-1,+1,-1,0.123,-0.234,1/8,   1 ]';
                QAr  = [3,4,2,5,6,7,8,9,10,1,11];
                ewid = 1.22;
            case 4,
                QAi  = [ 2, 4, 3, 6,3,         5,   7,  8 ]';
                QAj  = [ 1, 1, 2, 2,3,         3,   3,  3 ]';
                QAv  = [+1,-1,+1,-1,0.123,-0.234, 1/4,  1 ]';
                QAr  = [3,4,2,5,6,7,1,8];
                ewid = 0.84;
        end
        
        nQA     = max(QAi);
        mQA     = max(QAj);
        nc      = size(ndxs,2);
        new_n   = n + (nQA-3) * nc + 1;
        new_m   = m + mQA * nc;
        n_ndxs  = [ ndxs ; reshape( n + 1 : new_n - 1, nQA-3, nc  ) ];
        n_ndxs  = n_ndxs(QAr,:);
        if ~quiet,
            disp( sprintf( '   Approximation size: %d variables, %d equality constraints', new_n, new_m ) );
            disp( spacer );
        end

        % Stuff free variables into a lorentz cone to preserve warm start
        tfree = ones( 1, n );
        for k = 1 : length(cones),
            tfree(cones(k).indices) = 0;
        end
        tfree(xndxs) = 1;
        tfree = find(tfree);
        
        % Perform (x,y,z) ==> w transformation on A and C
        c (new_n,end) = 0;
        At(new_n,end) = 0;
        b (new_m,end) = 0;
        
        % Add new cone constraints
        lQA  = length(QAi);
        nc0  = 0:mQA:mQA*(nc-1);
        nc1  = ones(1,nc);
        At  = [ At, sparse( n_ndxs(QAi,:), ...
                QAj(:,nc1) + nc0(ones(lQA,1),:), ...
                QAv(:,nc1), new_n, mQA * nc ) ];
        
        endxs = n_ndxs(nQA-3,:) + (m+mQA-1+nc0) * new_n;
        fndxs = n_ndxs(3,:)     + (m+mQA-1+nc0) * new_n;
        
        ncone.type       = 'semidefinite';
        ncone.indices    = reshape(n_ndxs(1:nQA-2,:),3,(nQA-2)*nc/3);
        ncone(2).type    = 'nonnegative';
        ncone(2).indices = n_ndxs(nQA,:);
        ncone(3).type    = 'lorentz';
        ncone(3).indices = [ tfree(:) ; new_n ];
        cones(texp) = [];
        cones = [ cones, ncone ];

        amult = 1;
        rel_err = 1e-6;
        last_err = Inf;
        epow_i = 1 / epow;
        
        ninfeas  = 0;
        best_x   = [];
        best_px  = Inf;
        best_ox  = Inf;
        best_py  = Inf;
        best_oy  = -Inf;
        best_y   = [];
        solv_warn = false;
        use_init = false;
        if ~quiet,
            disp( ' Target     Conic    Solver' );
            disp( 'Precision   Error    Status' );
            disp( '---------------------------' );
        end
        best_x = NaN * ones(n,1);
        best_y = NaN * ones(m,1);
        dobj = -Inf;
        pobj = +Inf;
        stagnant = false;
        last_slow = false;
        last_u = false;
        dscale = blkdiag(speye(m,m),speye(new_m-m,new_m-m));
        XY0 = {};
        nprec = prec(3) * [ 1, 1, 1 ];
        for iter = 1 : 1000,
            x0e = x0 * epow_i;
            qq1 = exp( -x0e );
            qq2 = 1 - x0e;
            At(endxs) = - amult * qq1;
            At(fndxs) = amult * qq2;
            [ x, y, status, z ] = feval( sfunc, At, b, c, cones, true, nprec );
            if status(1) == 'F',
                tprec = Inf;
            else
                tprec = nprec(2+(status(3)=='a'));
            end
            tighten = false;
            if isnan(x(1)),
                x_valid = false;
                x_clean = false;
            else
                xxx = x(xndxs,:);
                yyy = x(yndxs,:);
                zzz = x(zndxs,:);
                x_valid = all( yyy >= 0 & zzz >= 0 );
                if x_valid,
                    yy2 = max( yyy, realmin );
                    nx0 = xxx ./ yy2;
                    nx1 = log( zzz ./ yy2 );
                    nxe = nx0 - nx1;
                    ttx = yyy == 0 | nxe <= 0;
                    x_clean = all( ttx );
                    pob = c' * x;
                else
                    x_clean = false;
                    tighten = true; 
                end
            end
            if isnan(y(1)),
                y_valid = false;
            else
                z = z + At * [ zeros(m,1) ; y(m+1:end) ];
                uuu = z( xndxs, : );
                vvv = z( yndxs, : );
                www = z( zndxs, : );
                y_valid = all( uuu <= 0 & www >= 0 );
                if y_valid,
                    uu2 = min( uuu, -realmin );
                    nx2 = 1 - vvv ./ uu2;
                    nx3 = - log( www ./ uu2 );
                    if all( uuu == 0 | nx3 <= nx2 ),
                        y_clean = true;
                        dob = b' * y;
                        nx0 = 0.5 * ( nx0 + nx2 );
                    else
                        y_valid = false;
                    end
                end
                if ~y_valid, 
                    tighten = true; 
                end
            end
            if x_valid & ~y_valid,
                nxe( ttx ) = 0;
                nx0( ttx ) = x0( ttx );
            end
            switch status,
                case { 'Infeasible', 'Inaccurate/Infeasible' },
                    x_clean = y_valid;
                    stagnant = true;
                    last_u = false;
                    pob = Inf;
                case { 'Unbounded', 'Inaccurate/Unbounded' },
                    y_valid = x_valid;
                    stagnant = true;
                    last_u = true;
                    dob = -Inf;
                case { 'Solved', 'Inaccurate/Solved' },
                    stagnant = ~last_u;
                    last_u = false;
            end
            if x_clean & ( tprec < best_px | ( tprec == best_px & pob < best_ox ) ),
                best_x  = x(1:n);
                best_px = tprec;
                best_ox = pob;
            end
            if y_valid & ( tprec < best_py | ( tprec == best_py & dob > best_oy ) ),
                best_y  = y(1:m);
                best_py = tprec;
                best_oy = dob;
            end
            if tighten,
                err = Inf;
            else
                err = max(0,max(nxe));
            end
            stagnant = stagnant & ~isnan( x(1) ) & ~isnan( y(1) ) & ( err > 0.9 * last_err );
            if ~quiet,
                if stagnant, stagc = 'S'; else stagc = ' '; end
                disp( sprintf( '%9.3e%c %9.3e  %s', nprec(2), stagc, err, status ) );
            end
            if err == 0,
                if nprec(1) == prec(1), break; end
                advance = true;
            elseif ~stagnant,
                advance = false;
            elseif tprec == nprec(2) & amult < 100,
                amult = amult * 100;
                At(:,m+1:end) = At(:,m+1:end) * 100;
            elseif nprec(1) == prec(1),
                break;
            else
                advance = true;
            end
            if advance,
                nprec(1:2) = prec(1:2);
                stagnant = false;
                err = Inf;
                At(:,m+1:end) = At(:,m+1:end) / amult;
                amult = 1;
            end
            last_err = err;
            x0 = min( max( nx0, x0 - epow ), x0 + epow );
            x0 = min( max( -maxw, x0 ), maxw );
        end
        if isnan( best_x(1) ), 
            status = 'Infeasible';
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
    else
        [ x, y, status ] = feval( sfunc, At, b, c, cones, quiet, prec );
    end
    if cvx___.profile, profile resume; end
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
    oval  = sgn * d;
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

gscale = 0;
if gobj,
    switch status,
        case 'Unbounded', gscale = 1; status = 'Solved';
        case 'Inaccurate/Unbounded',  gscale = 1; status = 'Inaccurate/Solved';
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
