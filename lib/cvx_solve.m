function [ status, result, bound, iters, tol ] = cvx_solve

global cvx___
try
    pstr = cvx___.problems(end);
catch
    cvx_throw( 'No CVX model is present.' );
end

quiet = pstr.quiet;
obj   = pstr.objective;
prec  = pstr.precision;
gobj  = abs( pstr.direction ) > 1;
solv  = pstr.solver;
shim  = cvx___.solvers.list( solv );
params = shim.eargs;
params.quiet = quiet;
params.precision = prec;
params.warmstart = cvx___.warmstart(2:end);
params.settings = pstr.settings;
ndual = ~isempty( pstr.duals );
nobj  = numel( obj );
if nobj > 1,
    cvx_throw( 'Your objective function is not a scalar.' );
end
[ At, cones, sgn, Q, P, exps, dualized ] = cvx_extract( shim.config, shim.name );
idual_error = false;
if ndual,
    ctype = { cones.type };
    if any( strcmp( ctype, 'integer' ) ) || any( strcmp( ctype, 'binary' ) ),
        idual_error = true;
        ndual = false;
    end
end

% Yes, the negative sign is here. This is new. I decided to make the 
% Lagrangian matrix in @cvxprob/eliminate fully consistent, which requires
% either the objective or the constraints to be negated.
c = - At( :, 1 );
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

zero_c = false ; % nnz( c ) == 0;
if zero_c,
    c = [ c ; 1 ]; %#ok
    At(end+1,:) = b * sqrt(mean(sum(At.^2))) / norm(b);
    cones(end+1) = struct( 'type', 'nonnegative', 'indices', n+1 );
    n = n + 1;
end

%
% Ferret out the degenerate and overdetermined problems
%

x     = NaN * ones(n,1);
y     = NaN * ones(m,1);
oval  = NaN;
bval  = NaN;
pval  = NaN;
dval  = NaN;
tprec = Inf;
estruc = [];
badsol = false;

iters = 0;
tt = ( b' ~= 0 ) & ~any( At, 1 );
infeas = any( tt );
if m > n && n > 0,
    
    %
    % Overdetermined problem
    %
    
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
    badsol = true;

elseif n ~= 0 && ~infeas && ( any( b ) || any( c ) ),
        
    %
    % Call solver
    %
    
    if isempty( cones ),
        texp = [];
    else
        texp = find( strcmp( { cones.type }, 'exponential' ) );
    end
    need_iter = ~isempty( texp ) && shim.config.dualize && ~isfield( shim.config, 'exponential' );
    cvx_setspath;
    if ~quiet,
        disp( ' ' );
        spacer = '-';
        if need_iter,
            disp( 'Successive approximation method to be employed.' );
        else
            sname = shim.name;
            if ~isempty( shim.version ), sname = [ sname, ' ', shim.version ]; end
            fprintf( 'Calling %s: %d variables, %d equality constraints\n', sname, n, m );
            spacer = spacer(:,ones(1,60));
        end
        if dualized,
            fprintf( '   For improved efficiency, %s is solving the dual problem.\n', shim.name );
        end
        if need_iter,
            fprintf( '   %s will be called several times to refine the solution.\n', shim.name );
            fprintf( '   Original size: %d variables, %d equality constraints\n', n, m );
            spacer = spacer(:,ones(1,65));
        else
            disp( spacer );
        end
    end
    if cvx___.profile, 
        pstat = profile('status');
        profon = isequal( pstat.ProfilerStatus, 'on' );
        profile off; 
    else
        profon = false;
    end
    tstart = tic;
    if need_iter,
        
        %
        % Cone:
        %     cl { (x,y,z) | y*exp(x/y) <= z, y > 0 }
        %   = cl { (x,y,z) | x <= -y*log(y/z), z > 0 }
        % Approximation: given a shift point x0,
        %    { (x,y,z) | y*exp(x0)*pos(1+(x/y-x0)/8)^8 <= z, y > 0 }
        % Transformed cone:
        %   3 lorentz cones, 1 slack
        %
        
        ndxs  = cat( 2, cones(texp).indices );
        nc    = size(ndxs,2);
        new_n = n + 7 * nc;
        new_m = m + 4 * nc;
        x0    = realmin * ones(nc,1);
        maxw  = log(realmax);
        ndx2  = [ ndxs ; reshape(n+1:n+7*nc,7,nc) ];
        ndx3  = 1 : n; 
        ndx3(ndxs) = [];
        params.quiet = true;
        
        epow = 8; epow_i = 0.125; g = 0.123; h = 0.345;
        Pr = [ 2, 3, 7,10, 2, 3, 2, 3]; Qr = [ 1, 5, 6, 2, 3, 5, 6, 2, 3, 8, 9, 4, 8, 9];
        Pc = [ 1, 1, 1, 1, 2, 2, 3, 3]; Qc = [ 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4];
        Pv = [ g,-g,+h,-1,-1, 1, 1, 1]; Qv = [-1, 1, 1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1];
        
        Pr = [vec(ndx2(Pr,:));ndx3'];
        Pc = [vec(ndx2(Pc,:));ndx3'];
        Pv = ones(nc,1) * Pv;
        Pw = ones(n-3*nc,1);
        QQ = sparse(ndx2(Qr,:),bsxfun(@plus,Qc',0:4:4*nc-1),Qv'*ones(1,nc),new_n,4*nc);
        amult = 1;
        b(new_m) = 0;
        
        if ~quiet,
            fprintf( '   %d exponentials add %d variables, %d equality constraints\n', nc, new_n - n, new_m - m );
            disp( spacer );
        end
        
        cones(texp) = [];
        cones = cvx_pushcone( cones, 'lorentz', reshape(ndx2(1:9,:),3,3*nc) );
        cones = cvx_pushcone( cones, 'nonnegative', ndx2(end,:) );
        
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
            
            ex0e = exp( - x0 * epow_i );
            ax = epow * ex0e;
            az = epow - x0;
            Pq = [az,-az,ax,Pv(:,4:end)]';
            PP = sparse(Pr,Pc,[Pq(:);Pw]);
            
            % Solve the approximation
            [ x, status, tprec, iters2, y ] = shim.solve( ...
                [ PP * At, bsxfun( @times, QQ, amult ) ], ...
                b, PP*c, cones, params );
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
                x = PP' * x;
                xxx = x(ndxs(1,:));
                yyy = max( realmin, x(ndxs(2,:)) );
                zzz = max( realmin, x(ndxs(3,:)) );
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
                z   = full(c - At * y(1:m));
                uuu = min( -realmin, z(ndxs(1,:),:) );
                vvv = z(ndxs(2,:),:);
                www = max( +realmin, z(ndxs(3,:),:) );
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
                cer = min( epow, max( max( abs( cxX .* acX ) ), max( abs( cxY .* acY ) ) ) );
            elseif x_valid,
                cxX = cxX .* acX;
                tlX = tlX .* acX;
                cer  = min( epow, max( max( abs( cxX .* acX ) ) ) );
            elseif y_valid,
                cxX = cxY .* acY;
                tlX = tlY .* acY;
                erX = erY;
                cer = min( epow, max( max( abs( cxY .* acY ) ) ) );
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
                if cer >= 0.9 * last_cer, 
                    stagc = 's'; 
                end
                if err >= 0.9 * last_err, 
                    stage = 's'; 
                end
            end
            if ~quiet,
                fprintf( '%3d/%3d | %9.3e%c %9.3e%c %9.3e | %s\n', nmov, nact, cer, stagc, err, stage, tol, status );
            end
            
            % Solution found or no more iterations
            % In perfect arithmetic, erY should be all zeros---because the
            % approximate dual should be feasible in the original, too. But
            % in imperfect arithmetic, it may not be. So, we're using that
            % error as a threshold to decide when the *primal* point is
            % sufficiently accurate, too.
            if found,
                if tprec(1) < best_prec,
                    best_x = x;
                    best_y = y;
                    best_prec = tprec(1);
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
                    if all( boost ), amult = amult * 10;
                    else amult = amult .* repmat(9*boost+1,[1,4]); end
                    amult = min( amult, 1e5 );
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
            badsol = true;
        elseif isnan( best_y(1) ), 
            status = 'Unbounded';
            badsol = true;
        else
            status = 'Solved';
        end
        if best_prec > prec(3),
            status = 'Failed';
            badsol = true;
        elseif best_prec > prec(2),
            status = [ 'Inaccurate/', status ];
        end
        x = best_x(1:n,:);
        y = best_y(1:m,:);
        c = c(1:n,:);
    elseif ndual || dualized,
        try
            [ x, status, tprec, iters, y ] = shim.solve( At, b, c, cones, params );
        catch estruc
            status = 'Error';
        end
    else
        try
            [ x, status, tprec, iters ] = shim.solve( At, b, c, cones, params );
        catch estruc
            status = 'Error';
        end
    end
    try
        cvx___.timers(4) = cvx___.timers(4) + ( tic - tstart );
    catch
        cvx___.timers(4) = double(cvx___.timers(4)) + ( double(tic) - double(tstart) );
    end
    if profon,
        profile resume; 
    end
    if ~cvx___.path.hold, 
        cvx_clearspath; 
    end
    if zero_c,
        q = x(end); %#ok
        x(end) = [];
        switch status,
        case { 'Solved', 'Inaccurate/Solved' },
            if q > prec(3),
                status = strrep( status, 'Solved', 'Infeasible' );
                oval = sgn * Inf;
                bval = oval;
                y = y / abs( b' * y );
                x(:) = NaN;
                dval = 0;
            else
                oval = 0;
                bval = 0;
                pval = 1;
                dval = 1;
            end
            badsol = true;
        otherwise,
            if ~isequal( status, 'Error' ), 
                status = Failed; 
            end
        end
    else
        switch status,
        case { 'Solved', 'Inaccurate/Solved', 'Suboptimal' },
            oval = sgn * ( c' * x + d' );
            if ndual || dualized,
                bval = sgn * ( b(1:m,:)' * y + d' );
            elseif length(tprec) > 1,
                bval = sgn * tprec(2) + d';
            else
                bval = sgn * -Inf;
            end
            pval = 1;
            dval = 1;
        case { 'Infeasible', 'Inaccurate/Infeasible' },
            badsol = true;
            oval = sgn * Inf;
            bval = oval;
            dval = 0;
        case { 'Unbounded', 'Inaccurate/Unbounded' },
            badsol = true;
            oval = -sgn * Inf;
            bval = oval;
            pval = 0;
        otherwise,
            badsol = true;
            bval = NaN;
            if ~isnan( x ), pval = 1; end
            if ~isnan( y ), dval = 1; end
        end
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
    b( ~tt ) = 0;
    y = - b / ( b' * b );
    oval = sgn * Inf;
    bval = oval;
    dval = 0;
    badsol = true;
    
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
    oval = sgn * d;
    bval = oval;
    pval = 1;
    dval = 1;
    badsol = true;
    
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

tol = tprec(1);
if badsol
    cvx___.warmstart(2:end) = [];
end

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
if ~isempty( exps ),
    tt = exps(:,3) == 1;
    cvx___.x( exps(tt,2) ) = min( 1e300, exp( cvx___.x( exps(tt,1) ) ) );
    tt = exps(:,3) ~= 1;
    cvx___.x( exps(tt,1) ) = log( cvx___.x( exps(tt,2) ) );
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
end
oval = full(oval);
bval = full(bval);
result = oval;
bound = bval;
if ~quiet,
    if length( oval ) == 1,
        fprintf( 'Optimal value (cvx_optval): %+g\n', oval );
    else
        fprintf( 'Optimal value (cvx_optval): (multiobjective)\n' );
    end
end

if isempty( estruc ) && idual_error,
    warning( 'CVX:IntegerDual', ...
[ 'Dual variables are not supported for problems involving integer variables.\n', ...
  'All dual variables were set to the value NaN.' ] );
end

if ~quiet,
    disp( ' '  );
end

if ~isempty( estruc ),
    cvx_throw( estruc );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
