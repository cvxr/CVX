function shim = cvx_mosek( shim )

% Copyright 2013 CVX Research, Inc.
% This source file a trade secret of CVX Research, Inc. and may not be 
% obtained or distributed without express written permission of the owner.

global cvx___
if ~isempty( shim.solve )
    return
end

fs = cvx___.fs;
mext = cvx___.mext;
mlen = length(mext);

fbase = 'mosekopt';
fname = [ fbase, '.', mext ];

old_dir = pwd;
cleanup = onCleanup(@()cd(old_dir));

is_new = isempty( shim.name );
if is_new
    shim.name = 'Mosek';
    shim.dualize = true;
    shim.version = 'unknown';
    fpaths = which( fbase, '-all' );
    oshim = shim;
    shim = [];
    if cvx___.cs, scmp = @strcmp; else scmp = @strcmpi; end
    for k = 1 : length(fpaths)
        fpath = fpaths{k};
        if ~scmp(fpath(end-mlen+1:end), mext), continue; end
        if any( scmp( fpath, fpaths(1:k-1) ) ), continue; end
        if ~exist( fpath, 'file' ), continue; end
        tshim = oshim;
        tshim.fullpath = fpath;
        shim = [ shim, tshim ]; %#ok
    end
    if isempty( shim )
        shim = oshim;
        if isempty( shim.error )
            shim.error = 'Could not find a MOSEK MEX file.';
        end
        return
    end
end

for k = 1 : length(shim)
    if ~isempty(shim(k).error)
        continue
    end
    fpath = shim(k).fullpath;
    fspos = strfind(fpath, fs);
    npath = fpath(1:fspos(end)-1);
    shim(k).location = npath;
    fpath = [ npath, fs, fname ];
    if ~exist( fpath, 'file' )
        shim(k).error = sprintf( 'The MOSEK MEX file expected at\n    %s\nseems to be missing.', fpath );
        continue
    end
    cd( npath );
    if is_new
        try
            otp = evalc(fbase);
        catch exc
            shim(k).error = sprintf( 'Unexpected MEX file failure:\n%s\n', exc.message );
            continue
        end
        otp = regexp( otp, 'MOSEK Version \S+', 'match' );
        if ~isempty(otp)
            shim(k).version = otp{1}(15:end);
            nversion = sum(sscanf(shim(k).version,'%d'));
        else
            nversion = 0;
        end
        if nversion < 7
            shims(k).error('The CVX/MOSEK interface requires MOSEK version 7 or later.');
            continue
        end
    end
    try
        [rr,res] = mosekopt('minimize echo(0)',struct('c',1,'a',sparse(1,1,1),'blc',0)); %#ok
        if res.rcode ~= 0
            emsg = { sprintf( 'Error code %d (%s).', res.rcode, res.rcodestr ) };
            shim(k).error = sprintf( '%s\n', emsg{:} );
            continue
        end
    catch errmsg
        shim(k).error = sprintf( 'Unexpected MEX file failure:\n%s\n', errmsg.message );
        continue
    end
    clear('mosekopt');
    shim(k).fullpath = fpath;
    shim(k).path = [ npath, cvx___.ps ];
    shim(k).check = @check;
    shim(k).solve = @solve;
    shim(k).eargs = { @mosekopt };
end

function found_bad = check( nonls, sdp, mfunc ) %#ok
found_bad = false;
if ~sdp
    for k = 1 : length( nonls )
        if any( strcmp( nonls(k).type, { 'semidefinite', 'hermitian-semidefinite' } ) ) && size(nonls(k).indices,1) > 4
            warning( 'CVX:Mosek:Semidefinite', 'This nonlinearity requires use of semidefinite cones which Mosek does not support.\n%s', ...
                'You will need to use a different solver for this model.' );
            found_bad = true;
        end
    end
end

function [ x, status, tol, iters, y, z ] = solve( At, b, c, nonls, quiet, prec, settings, mfunc )
zp = zeros(0,1);
[n,m] = size(At);
b = b(:); c = c(:);
prob  = struct( 'a', zp, 'blc', b, 'buc', b, 'blx', -Inf(n,1), 'bux', Inf(n,1), 'c', c, 'cones', {{}} );
prob.ints.sub = zp;
xscale = zp;
sdp_n = 0;
for k = 1 : length( nonls )
    nonl = nonls(k);
    ti = nonl.indices;
    nn = size( ti, 1 );
    nv = size( ti, 2 );
    need_sdp = false;
    switch ( nonl.type )
    case 'i_integer'
            
        prob.ints.sub = [ prob.ints.sub ; ti(:) ]; 
        
    case'i_binary'
        
        prob.ints.sub = [ prob.ints.sub ; ti(:) ]; 
        prob.blx( ti ) = 0; prob.bux( ti ) = 1;
        
    case 'i_semicontinuous'
        
        error( 'CVX:SolverIncompatible', 'MOSEK does not support semicontinous variables.' );
        
    case 'i_semiinteger'
        
        error( 'CVX:SolverIncompatible', 'MOSEK does not support semiinteger variables.' );
        
    case 'nonnegative'
        
        prob.blx( ti ) = 0;

    case 'lorentz'

        if nn == 1
            prob.blx( ti ) = 0;
        else
            ti = ti([end,1:end-1],:);
            for qq = 1 : nv
                prob.cones{end+1}.type = 'MSK_CT_QUAD';
                prob.cones{end}.sub = ti(:,qq)';
            end
        end

    case 'rotated-lorentz'

        if nn <= 2
            prob.blx( ti ) = 0;
        else
            ti = ti([end-1:end,1:end-2],:);
            for qq = 1 : nv
                prob.cones{end+1}.type = 'MSK_CT_RQUAD';
                prob.cones{end}.sub = ti(:,qq)';
            end
        end

    case 'semidefinite'

        if nn == 3
            ti = ti([1,3,2],:);
            xscale(ti(1:2,:)) = true;
            for qq = 1 : nv
                prob.cones{end+1}.type = 'MSK_CT_RQUAD';
                prob.cones{end}.sub  = ti(:,qq)';
            end
        else
            n2  = 0.5 * ( sqrt( 8 * nn + 1 ) - 1 );
            qq2 = ( 1 : nn )';
            qqq = qq2;
            col = ceil( ( n2 + 0.5 - 0.5 / n2 ) - sqrt( ( n2 + 0.5 )^2 - 2 * qqq ) );
            row = qqq - ( col - 1 ) .* ( 2 * n2 - col ) / 2;
            vv2 = 1 + ( row ~= col );
            vv1 = 1.0 ./ vv2;
            need_sdp = true;
        end

    case 'hermitian-semidefinite'

        if nn == 4
            ti = ti([1,4,2,3],:);
            xscale(ti(1:2,:)) = true;
            for qq = 1 : nv
                prob.cones{end+1}.type = 'MSK_CT_RQUAD';
                prob.cones{end}.sub  = ti(:,qq)';
            end
        else
            %   X >= 0 <==> exists [ Y1, Y2^T ; Y2, Y3 ] >= 0 s.t.
            %               Y1 + Y3 == real(X), Y2 - Y2^T == imag(X)
            % So: <C,X> = <CR,XR> + <CI,XI>
            %           = <CR,Y1+Y3> + <CI,Y2-Y2^T>
            %           = <CR,Y1> + <CI,Y2> + <CI^T,Y2^T> + <CR,Y3>
            %           = < [ CR, CI^T ; CI, CR ], [ Y1, Y2^T ; Y2, Y3 ] >
            n2   = sqrt( nn );
            qqq  = ( 1 : ( 2 * nn + n2 ) )';
            col  = ceil( ( 2 * n2 + 0.5 - 0.25 / n2 ) - sqrt( ( 2 * n2 + 0.5 )^2 - 2 * qqq ) );
            row  = qqq - ( col - 1 ) .* ( 4 * n2 - col ) / 2;
            row2 = rem( row - 1, n2 ) + 1;
            col2 = rem( col - 1, n2 ) + 1;
            dig  = row2 == col2;
            img  = row > n2 & col <= n2;
            neg  = img & row2 < col2;
            tmp  = col2(neg); col2(neg) = row2(neg); row2(neg) = tmp;
            qq2  = 2 * ( row2 - 1 ) + ( img | dig ) + ( 2 * n2 - 1 ) * ( col2 - 1 ) - col2 .* ( col2 - 1 );
            vv2  = ( 2 - dig ) .* ( 1 - 2 * neg );
            tt   = cumsum([n2+1,2*n2:-1:n2+2]);
            qqq(tt) = []; qq2(tt) = []; vv2(tt) = [];
            vv1  = 1.0 ./ vv2;
            vv2  = 0.5 * vv2;
            n2   = 2 * n2;
            nn   = n2 * ( n2 + 1 ) / 2;
            need_sdp = true;
        end

    case 'exponential'

        ti = ti([3,2,1],:);
        for qq = 1 : nv
            prob.cones{end+1}.type = 'MSK_CT_PEXP';
            prob.cones{end}.sub = ti(:,qq)';
        end

    otherwise

        error( 'Invalid cone type: %s', tt );

    end 
    if need_sdp

        if ~sdp_n
            prob.bardim = zp;
            prob.barc = struct( 'subj', zp, 'subk', zp, 'subl', zp, 'val', zp );
            prob.bara = struct( 'subi', zp, 'subj', zp, 'subk', zp, 'subl', zp, 'val', zp );
            sndxi = zp; sndxj = zp; sndxp = zp; sndxd = zp;
            nextmat = 1;
        end

        if nv > 1
            qqq = bsxfun( @plus, qqq, nn       * ( 0 : nv - 1 ) );
            qq2 = bsxfun( @plus, qq2, qq2(end) * ( 0 : nv - 1 ) );
            vv2 = repmat( vv2, [1,nv] );
            vv1 = repmat( vv1, [1,nv] );
        end
        qq2 = ti(qq2);
        F = sparse( qqq, qq2, vv1, qqq(end), length(c) );
        prob.bardim = [ prob.bardim ; n2(ones(1,nv),:) ];
        
        cc = F * c;
        if nnz( cc )
            [ rr, cc, vv ] = find( cc ); %#ok
            mat = floor((rr-1)/nn);
            rr  = rr - nn * mat;
            prob.barc.subj = [ prob.barc.subj ; mat + nextmat ];
            prob.barc.subk = [ prob.barc.subk ; row(rr) ];
            prob.barc.subl = [ prob.barc.subl ; col(rr) ];
            prob.barc.val  = [ prob.barc.val  ; vv ];
        end

        cc = F * At;
        if nnz( cc )
            [ rr, cc, vv ] = find( cc );
            mat = floor((rr-1)/nn);
            rr  = rr - nn * mat;
            prob.bara.subi = [ prob.bara.subi ; cc  ];
            prob.bara.subj = [ prob.bara.subj ; mat + nextmat ];
            prob.bara.subk = [ prob.bara.subk ; row(rr) ];
            prob.bara.subl = [ prob.bara.subl ; col(rr) ];
            prob.bara.val  = [ prob.bara.val  ; vv ];
        end
        
        sndxi = [ sndxi ; qq2(:) ]; %#ok
        sndxj = [ sndxj ; sdp_n + qqq(:) ]; %#ok
        sndxp = [ sndxp ; sign(vv2(:)) ]; %#ok
        sndxd = [ sndxd ; vv2(:) ]; %#ok
        sdp_n = sdp_n + nn * nv;
        nextmat = nextmat + nv;
        
    end
end
prob.a = At.';
if ~isempty( xscale )
    alpha = sqrt(2.0);
    xscale = xscale ~= 0;
    xscale(end+1:n) = false;
    prob.c(xscale) = prob.c(xscale) * alpha;
    prob.a(:,xscale) = prob.a(:,xscale) * alpha;
end
if sdp_n
    zndxs = 1 : n;
    zndxs(sndxi) = [];
    qndxs = zeros(1,n);
    qndxs(zndxs) = 1 : numel(zndxs);
    prob.c = prob.c(zndxs);
    prob.a = prob.a(:,zndxs);
    prob.blx = prob.blx(zndxs);
    prob.bux = prob.bux(zndxs);
    prob.ints.sub = qndxs(prob.ints.sub);
    for k = 1 : length( prob.cones )
        prob.cones{k}.sub = qndxs(prob.cones{k}.sub);
    end
end
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = prec(1);
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = prec(1);
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = prec(1);
param.MSK_DPAR_INTPNT_CO_TOL_INFEAS = prec(1);
param.MSK_DPAR_INTPNT_CO_TOL_MU_RED = prec(1);
param.MSK_DPAR_INTPNT_TOL_PFEAS = prec(1);
param.MSK_DPAR_INTPNT_TOL_DFEAS = prec(1);
param.MSK_DPAR_INTPNT_TOL_INFEAS = prec(1);
param.MSK_DPAR_INTPNT_TOL_REL_GAP = max(1e-14,prec(1));
command = sprintf( 'minimize info echo(%d)', 3*~quiet );
if isfield( settings, 'write' )
    wfile = settings.write;
    if ~ischar( wfile ) || ndims( wfile ) ~= 2 || size( wfile, 1 ) ~= 1, %#ok
        error( 'CVX:MosekError', 'write filename must be a string.' );
    end
    settings = rmfield( settings, 'write' );
    command = sprintf( '%s write(%s)', command, wfile );
end    
[ rr, res ] = cvx_run_solver( mfunc, command, prob, param, 'rr', 'res', settings, 3 ); %#ok
if isfield( res.sol, 'int' )
    sol = res.sol.int;
    has_dual = false;
elseif isfield( res.sol, 'bas' )
    sol = res.sol.bas;
    has_dual = true;
else
    sol = res.sol.itr;
    has_dual = true;
end
tol = prec(2);
if sdp_n
    zndxs = sparse( zndxs, 1:numel(zndxs), 1, n, numel(zndxs) );
    sndxp = sparse( sndxi, sndxj, sndxp, n, sdp_n );
    x = full( zndxs * sol.xx + sndxp * sol.barx );
else
    x = sol.xx;
end
if has_dual
    y = sol.slc - sol.suc;
    z = sol.slx - sol.sux;
    if isfield( sol, 'snx' )
        z = z + sol.snx; 
    end
    if sdp_n
        sndxd = sparse( sndxi, sndxj, sndxd, n, sdp_n );
        z = full( zndxs * z + sndxd * sol.bars );
    end
else
    y = NaN * ones(m,1);
    z = NaN * ones(n,1);
end
if ~isempty( xscale )
    x(xscale) = x(xscale) * alpha;
    z(xscale) = z(xscale) / alpha;
end
status = '';
lbound = [];
switch sol.solsta
    case { 'NEAR_PRIMAL_INFEASIBLE_CER', 'PRIMAL_INFEASIBLE_CER' }
        status = 'Infeasible';
        x(:) = NaN; scl = abs(b'*y); z = z / scl; y = y / scl;
        lbound = Inf;
    case { 'NEAR_DUAL_INFEASIBLE_CER', 'DUAL_INFEASIBLE_CER' }
        status = 'Unbounded';
        y(:) = NaN; z(:) = NaN; x = x / abs(c'*x);
        lbound = -Inf;
    case { 'OPTIMAL', 'NEAR_OPTIMAL', 'INTEGER_OPTIMAL', 'NEAR_INTEGER_OPTIMAL' }
        status = 'Solved';
        if has_dual
            lbound = b' * y;
        elseif res.info.MSK_IINF_MIO_NUM_RELAX > 0 && isfield( res.info, 'MSK_DINF_MIO_OBJ_BOUND' )
            lbound = res.info.MSK_DINF_MIO_OBJ_BOUND;
        end
    case 'PRIMAL_FEASIBLE'
        if res.info.MSK_IINF_MIO_NUM_RELAX > 0
            status = 'Suboptimal';
            if isfield( res.info, 'MSK_DINF_MIO_OBJ_BOUND' )
                lbound = res.info.MSK_DINF_MIO_OBJ_BOUND;
            end
        end
end
if isempty(status)
    switch sol.prosta
        case { 'ILL_POSED', 'PRIMAL_INFEASIBLE_OR_UNBOUNDED' }
            tol = Inf;
            x(:) = NaN; y(:) = NaN; z(:) = NaN; 
        case { 'PRIMAL_INFEASIBLE', 'NEAR_PRIMAL_INFEASIBLE' }
            status = 'Infeasible';
            x(:) = NaN; z = z / abs(b'*y); y = y / abs(b'*y);
            sol.solsta = sol.prosta;
        case { 'PRIMAL_AND_DUAL_INFEASIBLE' }
            status = 'Infeasible';
            x(:) = NaN; y(:) = NaN; z(:) = NaN; 
            sol.solsta = sol.prosta;
        case { 'DUAL_INFEASIBLE', 'NEAR_DUAL_INFEASIBLE' }
            status = 'Unbounded'; 
            y(:) = NaN; z(:) = NaN; x = x / abs(c'*x);
            sol.solsta = sol.prosta;
        otherwise
            if has_dual
                pobj = c' * x;
                dobj = b' * y;
                nrmc = norm( c );
                nrmb = norm( b );
                xc   = At' * x;
                zc   = At * y + z;
                if c' * x < 0
                    uerr = -nrmc * norm(xc) / max( nrmb, 1 ) / pobj;
                else
                    uerr = Inf;
                end
                if b' * y < 0
                    ferr = -nrmb * norm(zc) / max( nrmc, 1 ) / dobj;
                else
                    ferr = Inf;
                end
                perr = norm( xc - b ) /  ( 1 + nrmb );
                derr = norm( zc - c ) / ( 1 + norm( c ) );
                gerr = max( 0, pobj - dobj ) / max( 1, abs(pobj) );
                tol2 = max( [ perr, derr, gerr ] );
                tol  = min( [ tol2, ferr, uerr ] );
                if tol == tol2
                    status = 'Solved';
                    lbound = dobj;
                elseif tol == ferr
                    status = 'Infeasible';
                    x(:) = NaN; z = z / abs(dobj); y = y / abs(dobj);
                    lbound = Inf;
                else
                    status = 'Unbounded';
                    y(:) = NaN; z(:) = NaN; x = x / abs(pobj);
                    lbound = -Inf;
                end
            else
                warning( 'CVX:UnknownMosekStatus', 'Unknown MOSEK status: %s/%s', sol.prosta, sol.solsta );
                x(:) = NaN; y(:) = NaN; z(:) = NaN; 
                tol = Inf;
            end
    end
end
if tol > prec(3)
    status = 'Failed';
elseif strncmp( sol.solsta, 'NEAR_', 5 )
    tol = prec(3);
    status = [ 'Inaccurate/', status ];
elseif tol > prec(2)
    status = [ 'Inaccurate/', status ];
end
if ~isempty( lbound )
    tol(2) = lbound;
end
iters = 0;


