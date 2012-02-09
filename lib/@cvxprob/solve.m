function solve( prob )

global cvx___
p = index( prob );
pr = cvx___.problems(p);
nobj = numel(pr.objective);
if nobj > 1 && ~pr.separable,
    error( 'Your objective function is not a scalar.' );
end
quiet = cvx___.problems(p).quiet;
obj   = cvx___.problems(p).objective;
gobj  = cvx___.problems(p).geometric;
clear pr
[ At, cones, sgn, Q, P, dualized ] = eliminate( prob, true );

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
if m > n && n > 0,
    
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
        warning( [ 'CVX:', status ], estr );
    end
    pval = NaN;
    dval = NaN;

elseif n ~= 0 && ~infeas && ( any( b ) || any( c ) ),
        
    %
    % Call solver
    %
    
    prob = cvx___.problems( p );
    prec = prob.precision;
    if isempty( cones ),
        texp = [];
    else
        texp = find( strcmp( { cones.type }, 'exponential' ) );
    end
    need_iter = ~isempty( texp );
    solv = prob.solver;
    lsolv = lower(solv);
    sfunc  = [ 'cvx_solve_', lsolv ];
    cvx_setspath( solv );
    if ~quiet,
        disp( ' ' );
        spacer = '-';
        if need_iter,
            disp( 'Successive approximation method to be employed.' );
        else
            fprintf( 1, 'Calling %s: %d variables, %d equality constraints\n', solv, n, m );
            spacer = spacer(:,ones(1,60));
        end
        if dualized,
            fprintf( 1, '   For improved efficiency, %s is solving the dual problem.\n', solv );
        end
        if need_iter,
            fprintf( 1, '   %s will be called several times to refine the solution.\n', solv );
            fprintf( 1, '   Original size: %d variables, %d equality constraints\n', n, m );
            spacer = spacer(:,ones(1,65));
        else
            disp( spacer );
        end
    end
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
                % ewid = 1.75;
            case 8,
                QAi  = [ 2, 4, 3, 6, 5, 7, 6, 9,3,         8,  10, 11 ]';
                QAj  = [ 1, 1, 2, 2, 3, 3, 4, 4,5,         5,   5,  5 ]';
                QAv  = [+1,-1,+1,-1,+1,-1,+1,-1,1.001,-1.002, 1/8,  1 ]';
                QAr  = [3,4,2,5,6,7,8,9,10,1,11];
                % ewid = 1.22;
            case 4,
                QAi  = [ 2, 4, 3, 6,3,         5,   7,  8 ]';
                QAj  = [ 1, 1, 2, 2,3,         3,   3,  3 ]';
                QAv  = [+1,-1,+1,-1,1.001,-1.002, 1/4,  1 ]';
                QAr  = [3,4,2,5,6,7,1,8];
                % ewid = 0.84;
        end
        
        nQA     = max(QAi);
        mQA     = max(QAj);
        nc      = size(ndxs,2);
        new_n   = n + (nQA-3) * nc; % + 1;
        new_m   = m + mQA * nc;
        n_ndxs  = [ ndxs ; reshape( n + 1 : new_n, nQA-3, nc  ) ];
        n_ndxs  = n_ndxs(QAr,:);
        if ~quiet,
            fprintf( 1, '   %d exponentials add %d variables, %d equality constraints\n', nc, new_n - n, new_m - m );
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
        cones(texp) = [];
        cones = [ cones, ncone ];
        
        arow = size(Anew,2) / nc;
        orow = ones(arow,1);
        amult = ones(1,nc);
        epow_i = 1 / epow;

        oprec = prec;
        best_x = NaN * ones(n,1);
        best_y = NaN * ones(m,1);
        best_prec = Inf;
        if ~quiet, 
            disp( ' Cones  |             Errors              |' );
            disp( 'Mov/Act | Centering  Exp cone   Poly cone | Status' );
            disp( '--------+---------------------------------+---------' );
        end
        failed = 0;
        attempts = 0;
        last_err = Inf;
        last_cer = Inf;
        last_solved = 0;
        max_eiters = 25;
        for iter = 1 : max_eiters,
            % Insert the current centerpoints into the constraints
            x0e = x0 * epow_i;
            ex0e = exp( -x0e );
            Anew(endxs) = - ex0e;
            Anew(fndxs) = 1 - x0e; 
            Anew2 = Anew * diag(sparse(vec(amult(orow,:))));
            
            % Solve the approximation
            [ x, y, status, tprec, iters2, z ] = feval( sfunc, [ At, Anew2 ], [ b ; bnew ], c, cones, true, prec );
            iters = iters + iters2;
            x_valid = ~any(isnan(x));
            y_valid = ~any(isnan(y));
           
            % The approximate primal cone is a strict subset of the exact
            % primal cone. A point that is feasible for the approximate 
            % model is guaranteed to be feasible for the exact model only
            % if x/y == x0. Furthermore, the larger |x/y-x0| is, the weaker
            % the approximation. So, our goal in these iterations is to
            % minimize |x/y - x0|. The hope is that we can reduce this gap
            % to the point that the deviation from exact feasibility is
            % within our desired numerical tolerance.
            % Exact:  y .* exp( x ./ y ) <= z
            % Approx: exp(x0) .* y .* max(0,1+(x./y-x0)/p).^p <= z
            if x_valid,
                xxx = x( xndxs, : );
                yyy = max( realmin, x( yndxs, : ) );
                zzz = max( realmin, x( zndxs, : ) );
                nmX = sqrt( xxx .^ 2 + yyy .^ 2 + zzz .^ 2 );
                xxy = xxx ./ yyy - x0;
                zzy = zzz ./ yyy;
                xxz = log( zzy ) - x0;
                xxc = epow * ( max( 0, zzy .^ epow_i .* ex0e ) - 1 );
                tlX = max( 0, xxy - xxc );
                erX = max( 0, xxy - xxz );
                acX = erX ~= 0;
                cxX = xxy + 0.5 * ( xxz - xxy ) .* acX;
                ttt = yyy == realmin;
                if any( ttt ),
                    xxy( ttt ) = -2 * ( sign( xxx( ttt ) ) * realmax );
                    cxX( ttt ) = max( 1 - epow, min( xxy( ttt ), epow - 1 ) );
                    tlX( ttt ) = 0;
                    erX( ttt ) = 0;
                end
            end
            
            % The exact dual cone is a strict subset of the approximate
            % dual cone. Therefore any point that is dual feasible in the
            % approximate model is also dual feasible in the exact model.
            % The further x/y is from x0, the farther away such a point
            % will be from the boundary of the exact dual cone.
            % Exact:  -u.*exp(v/u-1)<=w
            % Approx: -exp(-x0).*u.*(1-(v./u-1+x0)/(p-1)).^(1-p)<=w
            if y_valid,
                z = z + Anew2 * y(m+1:end);
                uuu = min( -realmin, z( xndxs, : ) );
                vvv = z( yndxs, : );
                www = max( +realmin, z( zndxs, : ) );
                nmY = sqrt( uuu .^ 2 + vvv .^ 2 + www .^ 2 );
                wwu = - www ./ uuu;
                xxu = 1 - vvv ./ uuu - x0;
                xxw = - log( wwu ) - x0;
                xxd = ((exp(x0).*wwu).^(1/(1-epow))-1)*(epow-1);
                tlY = max( 0, xxd - xxu );
                erY = max( 0, xxw - xxu );
                acY = erY ~= 0;
                cxY = xxu + 0.5 * ( xxw - xxu ) .* acY;
                ttt = uuu == -realmin;
                if any( ttt ),
                    cxY( ttt ) = max( 1 - epow, min( xxu( ttt ), epow - 1 ) );
                    tlY( ttt ) = 0;
                    erY( ttt ) = 0;
                end
            end
            
            if x_valid && y_valid,
                cxX = ( nmX .* cxX + nmY .* cxY ) ./ ( nmX + nmY );
                kkt = ( xxx .* uuu + yyy .* vvv + zzz .* www ) ./ ( nmX .* nmY ) < 1e-4;
                kkX = ( nmX > 1e-3 * nmY ) | kkt;
                kkY = ( nmY > 1e-3 * nmX ) | kkt;
                erX = erX .* kkX;
                tlX = tlX .* kkX;
                acX = acX .* kkX;
                acY = acY .* kkY;
                cer  = min( epow, max( max( abs( cxX .* acX ) ), max( abs( cxY .* acY ) ) ) );
            elseif x_valid,
                cxX = cxX .* acX;
                tlX = tlX .* acX;
                cer  = min( epow, max( max( abs( cxX .* acX ) ) ) );
            elseif y_valid,
                cxX = cxY .* acY;
                tlY = tlY .* acY;
                erX = erY;
                cer  = min( epow, max( max( abs( cxY .* acY ) ) ) );
            end
            if x_valid || y_valid,
                err  = max( erX );
                tol  = max( tlX );
                nmov = nnz( erX > max( prec(2), 1.5 * tlX ) );
                nact = nnz( erX );
            else
                err = 0; tol = 0; cer = 0;
                nmov = 0; nact = 0;
            end
            solved = x_valid * 2 + y_valid;
            found = nmov == 0 && solved;
                
            
            % Check for stagnation
            stagc = ' '; stage = ' ';
            if ~found && last_solved == solved && last_act == nact,
                if cer >= 0.9 * last_cer, stagc = 's'; end
                if err >= 0.9 * last_err, stage = 's'; end
            end
            if ~quiet,
                fprintf( 1, '%3d/%3d | %9.3e%c %9.3e%c %9.3e | %s\n', nmov, nact, cer, stagc, err, stage, tol, status );
            end
            
            % Solution found or no more iterations
            % In perfect arithmetic, erY should be all zeros---because the
            % approximate dual should be feasible in the original, too. But
            % in imperfect arithmetic, it may not be. So, we're using that
            % error as a threshold to decide when the *primal* point is
            % sufficiently accurate, too.
            if found,
                if tprec < best_prec,
                    best_x = x;
                    best_y = y;
                    best_prec = tprec;
                end
                if best_prec <= prec(1) || attempts == 2,
                    break;
                end
                attempts = attempts + 1;
            end
            if status(1) == 'F',
                failed = failed + 0.5 * ( 1 + ~x_valid );
                if failed >= 3, break; end
                if ~x_valid,
                    prec(3) = prec(3) * 10;
                    continue;
                end
            else
                prec(3) = oprec(3);
                failed = 0;
            end
            
            % Stagnation?
            if stagc == 's' || stage == 's',
                if all( amult == 1e5 ), break; end
                amult = min( amult * 10, 1e5 ); 
            elseif ~failed,
                boost = ~cxX & erX;
                if any( boost ),
                    amult(boost) = min( amult(boost) * 10, 1e5 );
                end
            end
            
            % Shift centerpoint
            last_solved = x_valid * 2 + y_valid;
            last_cer = cer;
            last_err = err;
            last_act = nact;
            if last_solved,
                x0 = max( min( x0 + max( min( epow, cxX ), -epow ), maxw ), -maxw );
            end
            
        end
        if isnan( best_x(1) ), 
            status = 'Infeasible';
        elseif isnan( best_y(1) ), 
            status = 'Unbounded';
        else
            status = 'Solved';
        end
        if best_prec > prec(3),
            status = 'Failed';
        elseif best_prec > prec(2),
            status = [ 'Inaccurate/', status ];
        end
        x = best_x(1:n,:);
        y = best_y(1:m,:);
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
    tprec = 0;
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
    tprec = 0;
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

trick = false;
if gobj,
    switch status,
        case 'Unbounded', 
            status = 'Solved';
            trick = true;
        case 'Inaccurate/Unbounded', 
            status = 'Inaccurate/Solved';
            trick = true;
    end
end

if ~quiet,
    fprintf( 1, 'Status: %s\n', status );
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
    if trick, y = P(:,1) + realmax * sign(y); end
    cvx___.x = y;
    cvx___.y = x(2:end);
else
    if trick, x = Q(:,1) + realmax * sign(x); end
    cvx___.x = x;
    cvx___.y = y(2:end);
end
if nnz( cvx___.exponential ),
    esrc = find( cvx___.exponential );
    edst = cvx___.exponential( esrc );
    cvx___.x( edst ) = min( 1e300, exp( cvx___.x( esrc ) ) );
end

%
% Compute the objective
%

if ~isempty( obj ),
    if isinf( oval ) || isnan( oval ),
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
        fprintf( 1, 'Optimal value (cvx_optval): %+g\n', oval );
    else
        fprintf( 1, 'Optimal value (cvx_optval): (multiobjective)\n' );
    end
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
