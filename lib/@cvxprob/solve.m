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

iters = 0;
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
                QAv  = [+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,1.001,-1.002,1/16,  1 ]';
                QAr  = [3,4,2,5,6,7,8,9,10,11,12,13,1,14];
                ewid = 1.75;
            case 8,
                QAi  = [ 2, 4, 3, 6, 5, 7, 6, 9,3,         8,  10, 11 ]';
                QAj  = [ 1, 1, 2, 2, 3, 3, 4, 4,5,         5,   5,  5 ]';
                QAv  = [+1,-1,+1,-1,+1,-1,+1,-1,1.001,-1.002, 1/8,  1 ]';
                QAr  = [3,4,2,5,6,7,8,9,10,1,11];
                ewid = 1.22;
            case 4,
                QAi  = [ 2, 4, 3, 6,3,         5,   7,  8 ]';
                QAj  = [ 1, 1, 2, 2,3,         3,   3,  3 ]';
                QAv  = [+1,-1,+1,-1,1.001,-1.002, 1/4,  1 ]';
                QAr  = [3,4,2,5,6,7,1,8];
                ewid = 0.84;
        end
        
        nQA     = max(QAi);
        mQA     = max(QAj);
        nc      = size(ndxs,2);
        new_n   = n + (nQA-3) * nc; % + 1;
        new_m   = m + mQA * nc;
        n_ndxs  = [ ndxs ; reshape( n + 1 : new_n, nQA-3, nc  ) ];
        n_ndxs  = n_ndxs(QAr,:);
        if ~quiet,
            disp( sprintf( '   Approximation size: %d variables, %d equality constraints', new_n, new_m ) );
            disp( spacer );
        end

        % Stuff free variables into a lorentz cone to preserve warm start
        % tfree = ones( 1, n );
        % for k = 1 : length(cones),
        %     tfree(cones(k).indices) = 0;
        % end
        % tfree(xndxs) = 1;
        % tfree = find(tfree);
        
        % Perform (x,y,z) ==> w transformation on A and C
        c (new_n,end) = 0;
        At(new_n,end) = 0;
        
        % Add new cone constraints
        lQA   = length(QAi);
        nc0   = 0:mQA:mQA*(nc-1);
        nc1   = ones(1,nc);
        Anew  = sparse( n_ndxs(QAi,:), ...
                 QAj(:,nc1) + nc0(ones(lQA,1),:), ...
                 QAv(:,nc1), new_n, mQA * nc );
        bnew  = zeros( new_m - m, 1 );
        
        endxs = n_ndxs(nQA-3,:) + (mQA-1+nc0) * new_n;
        fndxs = n_ndxs(3,:)     + (mQA-1+nc0) * new_n;
        
        ncone.type       = 'semidefinite';
        ncone.indices    = reshape(n_ndxs(1:nQA-2,:),3,(nQA-2)*nc/3);
        ncone(2).type    = 'nonnegative';
        ncone(2).indices = n_ndxs(nQA,:);
%       ncone(3).type    = 'lorentz';
%       ncone(3).indices = [ tfree(:) ; new_n ];
        cones(texp) = [];
        cones = [ cones, ncone ];

        amult = 1; bmult = 1; 
        epow_i = 1 / epow;
        
        best_x   = NaN * ones(n,1);
        best_px  = Inf;
        best_ox  = Inf;
        best_py  = Inf;
        best_oy  = -Inf;
        best_y   = NaN * ones(m,1);
        if ~quiet,
            disp( '    Precision         Conic     Solver' );
            disp( ' Target   Achieved    Error     Status' );
            disp( '-----------------------------------------' );
        end
        dobj = -Inf;
        pobj = +Inf;
        failed = 0;
        stagnant = 0;
        attempts = 0;
        success_once = false;
        last_nxe = 0;
        last_err = Inf;
        last_centered = 0;
        last_solved = false;
        nprec = prec(3) * [ 1, 1, 1 ];
        for iter = 1 : 1000,
            
            % Insert the current centerpoints into the constraints
            x0e = x0 * epow_i;
            Anew(endxs) = - exp( -x0e );
            Anew(fndxs) = 1 - x0e;
            
            % Solve the approximation
            [ x, y, status, tprec, iters2, z ] = feval( sfunc, bmult * [ At, amult * Anew ], bmult * [ b ; bnew ], c, cones, true, nprec );
            iters = iters + iters2;
            
            % Exit if we have nothing to work with
            x_isnan = any(isnan(x));
            y_isnan = any(isnan(y));
            if x_isnan && y_isnan || status(1) == 'F' && success_once,
                break;
            end

            % Tweak the scaling
            if ~x_isnan,
                norm_1 = norm( Anew' * x );
                norm_2 = norm( b * ~y_isnan - At' * x );
            end
            
            % The dual point should be feasible at all times
            err = Inf;
            y_valid = false;
            if ~y_isnan,
                yold = y(1:m);
                z = z + amult * Anew * y(m+1:end);
                uuu = z( xndxs, : );
                vvv = z( yndxs, : );
                www = z( zndxs, : );
                y_valid = all( uuu <= 0 & www >= 0 );
                if y_valid,
                    uu2 = min( uuu, -realmin );
                    nx2 = min( max( 1 - vvv ./ uu2, -maxw ), maxw );
                    nx3 = min( max( - log( www ./ uu2 ), -maxw ), maxw );
                    y_valid = all( uuu == 0 | nx3 <= nx2 );
                    if y_valid, 
                        if any(strfind(status,'Infeasible')), 
                            err = 0;
                            dob = +Inf;
                        else
                            dob = b' * yold; 
                        end
                        if tprec < best_py || ( tprec == best_py && dob > best_oy ),
                            best_y  = yold;
                            best_py = tprec;
                            best_oy = dob;
                        end
                    end
                end
            end
            
            % The primal point is only feasible when we're close enough
            x_valid = false;
            if ~x_isnan,
                xxx = x( xndxs, : );
                yyy = x( yndxs, : );
                zzz = x( zndxs, : );
                x_valid = all( yyy >= 0 & zzz >= 0 );
                if x_valid,
                    yy2 = max( yyy, realmin );
                    nx0 = xxx ./ yy2;
                    nx1 = log( zzz ./ yy2 );
                    nxe = max( 0, ( nx0 - nx1 ) .* ( yyy > 0 ) );
                    n_viol = nnz( nxe );
                    err = max( nxe );
                    x_valid = true;
                    if err == 0,
                        if any(strfind(status,'Unbounded')), pob = -Inf;
                        else pob = c' * x; end
                        if tprec < best_px || ( tprec == best_px && pob < best_ox ),
                            best_x  = x(1:n);
                            best_px = tprec;
                            best_ox = pob;
                        end
                    end
                end
            end
            
            % Determine new centerpoint
            is_stagnant = false;
            is_failed = false;
            advance = false;
            if x_valid,
                if y_valid, success_once = true; end
                tt1 = nx0 >= x0 & nx1 <= x0 & nxe;
                tt2 = nx1 > x0 & nxe;
                tt3 = nx0 < x0 & nxe;
                n_centered = nnz( tt1 );
                n_lost = nnz( nxe & ~last_nxe );
                last_nxe = nxe;
                if n_centered == n_viol,
                    nx0 = 0.5 * ( nx0 + nx1 );
                    is_stagnant = n_viol ~= 0 && n_lost == 0;
                else
                    tt2 = nx1 >= x0 & nxe;
                    nx0(tt1) = 0.5 * ( nx0(tt1) + nx1(tt1) );
                    nx0(tt2) = nx1(tt2);
                    if y_valid,
                        is_stagnant = last_solved && err > 0.5 * last_err && n_lost == 0;
                        last_solved = true;
                        last_err = err;
                    else
                        last_solved = false;
                    end
                end
            elseif y_valid,
                last_solved = false;
                n_centered = 0;
                nx0 = nx2;
            else
                is_failed = true;
            end
            
            % Stagnation?
            if ~is_stagnant,
                stagnant = 0;
            elseif stagnant == 2,
                advance = true;
            else
                stagnant = stagnant + 1;
                amult = amult * 100;
            end
            
            % Several invalid points in a row?
            if ~is_failed,
                failed = 0;
            elseif failed == 2,
                if x_isnan && y_isnan, break;
                else advance = true; end
            else
                failed = failed + 1;
            end
            
            % Print status
            if ~quiet,
                if is_stagnant, stagc = 'S'; else stagc = ' '; end
                disp( sprintf( '%9.3e %9.3e %9.3e%c  %s', nprec(2), tprec, err, stagc, status ) );
            end

            % Time to boost precision?
            if err == 0,
                advance = tprec <= nprec(2) || attempts == 1;
                attempts = attempts + 1;
            end
            if advance,
                if min(tprec,nprec(2)) <= prec(2), break; end
                nprec(1:2) = prec(1:2);
                last_solved = false;
                last_err = Inf;
                attempts = 0;
                stagnant = 0;
                failed = 0;
            end
            
            % Shift centerpoint
            if x_valid | y_valid,
                x0 = min( max( nx0, x0 - epow ), x0 + epow );
                x0 = min( max( -maxw, x0 ), maxw );
            end
        end
        if isnan( best_x(1) ), 
            status = 'Infeasible';
        elseif isnan( best_y(1) ), 
            status = 'Unbounded';
        else
            status = 'Solved';
        end
        tprec = max( best_px, best_py );
        if tprec > prec(3),
            status = 'Failed';
        elseif tprec > prec(2),
            status = [ 'Inaccurate/', status ];
        end
        x = best_x;
        y = best_y;
        c = c(1:n,:);
    else
        [ x, y, status, tprec, iters ] = feval( sfunc, At, b, c, cones, quiet, prec );
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
        if isnan( x ), pval = NaN; else pval = 1; end
        if isnan( y ), dval = NaN; else dval = 1; end
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
cvx___.problems( p ).iters = iters;
cvx___.problems( p ).tol = tprec;

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
