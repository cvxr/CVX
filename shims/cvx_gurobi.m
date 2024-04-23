function shim = cvx_gurobi( shim )

% Copyright 2018 CVX Research, Inc.
% This source file a trade secret of CVX Research, Inc. and may not be 
% obtained or distributed without express written permission of the owner.

global cvx___
if ~isempty( shim.solve )
    return
end

fs = cvx___.fs;
mext = cvx___.mext;
mlen = length(mext);

fbase = 'gurobi';
fname = [ 'gurobi.', mext ];

old_dir = pwd;
cleanup = onCleanup(@()cd(old_dir));

is_new = isempty( shim.name );
if is_new
    shim.name = 'Gurobi';
    shim.dualize = false;
    shim.version = 'unknown';
    fpaths = which( fbase, '-all' );
    switch mext
        case 'mexmaca64', d1 = '/Library/gurobi*/*/matlab';
        case 'mexmaci64', d1 = '/Library/gurobi*/*/matlab';
        case 'mexa64',    d1 = '/opt/gurobi*/*/matlab';
        case 'mexw64',    d1 = 'C:\gurobi*\*\matlab';
        otherwise,        d1 = '';
    end
    temp = dir( [ d1, fs, fname ] );
    for k = 1:length(temp)
        fpaths{end+1} = strcat( temp(k).folder, fs, temp(k).name ); %#ok
    end
    temp = dir( [ getenv( 'GUROBI_HOME' ), fs, 'matlab', fs, fname ] );
    for k = 1:length(temp)
        fpaths{end+1} = strcat( temp(k).folder, fs, temp(k).name ); %#ok
    end
    oshim = shim;
    shim = [];
    if cvx___.cs, scmp = @strcmp; else, scmp = @strcmpi; end
    for k = 1 : length(fpaths)
        fpath = fpaths{k};
        if any( scmp( fpath, fpaths(1:k-1) ) ), continue; end
        if ~exist( fpath, 'file' ), continue; end
        tshim = oshim;
        tshim.fullpath = fpath;
        shim = [ shim, tshim ]; 
    end
    if isempty( shim )
        shim = oshim;
        if isempty( shim.error )
            shim.error = 'Could not find a Gurobi MEX file.';
        end
        return
    end
end

prob = struct( 'Obj', 1, 'A', sparse(1,1,1), 'Sense', '>', 'RHS', 0 );
params = struct( 'OutputFlag', 0 );

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
        shim(k).error = sprintf( 'The Gurobi MEX file expected at\n    %s\nseems to be missing.', fpath );
        continue
    end
    cd( npath );
    try
        res = gurobi( prob, params );
    catch errmsg
        if any( strfind( errmsg.message, 'No valid Gurobi license was found.') ) || ...
           any( strfind( errmsg.message, 'No unlock license found.' ) )
            emsg = { 'No valid Gurobi license was found.' };
            shim(k).error = sprintf( '%s\n', emsg{:} );
        else
            shim(k).error = errmsg.message;
        end
        continue
    end
    vi = res.versioninfo;
    if vi.major < 9
        shim(k).error = 'CVX requires Gurobi 9.0 or later.';
        continue
    end
    shim(k).version = sprintf( '%d.%d%d', vi.major, vi.minor, vi.technical );
    shim(k).fullpath = fpath;
    shim(k).check = @check;
    shim(k).solve = @solve;
    shim(k).path = [ npath, cvx___.ps ];
    shim(k).eargs = {};
end

% GUROBI_CHECK
%
% We don't actually use this yet. The intention is to provide a gentle
% warning to the user if he selects the Gurobi solver while a model is
% being built; AND if the nonlinearities already present in the model are
% incompatible. We can also use it to provide warnings as constraints are
% entered. For now, however, we simply wait until the GUROBI_SOLVE sees 
% the problem and exits with a fatal error.

function found_bad = check( nonls )
found_bad = false;
for k = 1 : length( nonls )
    if any( strcmp( nonls(k).type, { 'semidefinite', 'hermitian-semidefinite' } ) ) && size(nonls(k).indices,1) > 4
        warning( 'CVX:SolverIncompatible', ...
            [ 'Gurobi does not support semidefinite cones larger than 2x2.\n', ...
              'You will need to use a different solver for this model.' ] );
        found_bad = true;
        break;
    end
end

% GUROBI_SOLVE
%
% This routine accepts the problem to solve in internal CVX form and
% performs the conversions necessary for Gurobi to solve it.

function [ x, status, tol, iters, y, z ] = solve( At, b, c, nonls, quiet, prec, settings )
need_y = nargout > 4;
need_z = nargout > 5;

n = numel(c);
m = numel(b);
prob       = [];
prob.obj   = full(c);
prob.A     = At';
prob.rhs   = full(b);
prob.sense = '=';
prob.lb    = -Inf * ones(n,1);
prob.ub    = -prob.lb;
prob.quadcon = struct( 'Qrow', {}, 'Qcol', {}, 'Qval', {}, 'q', {}, 'rhs', {} );
qz = zeros(n, 1);
prob.vtype = 'C';
prob.vtype = prob.vtype(1,ones(1,n));
is_int = false;
for k = 1 : length( nonls )
    nonl = nonls(k);
    tt = nonl.type;
    ti = nonl.indices;
    [ ni, mi ] = size( ti );
    qv = [];
    switch tt
        case 'i_integer'
            prob.vtype( ti ) = 'I';
            is_int = true;
        case 'i_binary'
            prob.vtype( ti ) = 'B';
            prob.lb( ti ) = 0;
            prob.ub( ti ) = 1;
            is_int = true;
        case 'i_semicontinuous'
            prob.vtype( ti ) = 'S';
            is_int = true;
        case 'i_semiinteger'
            prob.vtype( ti ) = 'N';
            is_int = true;
        case 'nonnegative'
            prob.lb( ti ) = 0;
        case 'lorentz'
            prob.lb( ti(end, :) ) = 0;
            if ni > 1
                t1 = ti;
                t2 = ti;
                qv = [ones(ni - 1, 1); -1];
            end
        case 'semidefinite'
            prob.lb( ti([1, end], :) ) = 0;
            if ni > 3
                error( 'CVX:SolverIncompatible', 'Gurobi does not support semidefinite cones larger than 2x2.\nYou must use another solver for this problem.' );
            elseif ni == 3
                t1 = ti([1,2], :);
                t2 = ti([3,2], :);
                qv = [-1, 1];
            end
        case 'hermitian-semidefinite'
            prob.lb( ti([1, end], :) ) = 0;
            if ni > 4
                error( 'CVX:SolverIncompatible', 'Gurobi does not support semidefinite cones larger than 2x2.\nYou must use another solver for this problem.' );
            elseif ni == 4
                t1 = ti([1,2,3], :);
                t2 = ti([4,2,3], :);
                qv = [-1, 1, 1];
            end
        case 'exponential'
            error( 'CVX:SolverIncompatible', 'Gurobi does not support the exponential cone.\nYou must use another solver for this problem.' );
        otherwise
            error( 'Invalid cone type: %s', tt );
    end
    if ~isempty(qv)
        for qq = 1 : mi
            prob.quadcon(end+1) = struct('Qrow', t1(:,qq), 'Qcol', t2(:,qq), 'Qval', qv, 'q', qz, 'rhs', 0);
        end
    end
end
prec(1) = prec(2);
params.OutputFlag = double(~quiet);
params.InfUnbdInfo = 1;
params.QCPDual = double(need_y);
params.BarConvTol = prec(1);
params.BarQCPConvTol = prec(1);
params.FeasibilityTol = max([1e-9,prec(1)]);
params.OptimalityTol = max([1e-9,prec(1)]);
res = cvx_run_solver( @gurobi, prob, params, 'res', settings, 3 );
tol = prec(2);
x = []; y = []; z = [];
lbound = [];
switch res.status
    case { 'NUMERIC', 'INF_OR_UNBD' }
        tol = Inf;
    case 'INFEASIBLE'
        status = 'Infeasible';
        lbound = -Inf;
        if isfield( res, 'farkasdual' )
            y = - res.farkasdual / abs( prob.rhs' * res.farkasdual );
            if need_z, z = prob.A' * y; end
        elseif need_y
            tol = Inf;
        end
    case 'UNBOUNDED'
        status = 'Unbounded';
        lbound = Inf;
        if isfield( res, 'unbdray' )
            x = res.unbdray / abs( prob.obj' * res.unbdray );
        else
            tol = Inf;
        end
    case { 'OPTIMAL', 'SUBOPTIMAL', 'INTERRUPTED', 'TIME_LIMIT' }
        status = 'Solved';
        if isfield( res, 'x' )
            x = res.x;
        else
            tol = Inf;
        end
        if isfield( res, 'pi' )
            y = res.pi;
            if need_z, z = prob.obj - prob.A' * y; end
            lbound = prob.rhs' * y;
        elseif need_y
            tol = Inf;
        end
        if isfield( res, 'objbound' )
            lbound = res.objbound;
        end
        if tol == Inf
            status = 'Failed';
        elseif res.status(1) == 'S' && ~is_int
            status = 'Inaccurate/Solved';
        elseif res.status(1) ~= 'O'
            status = 'Suboptimal';
        end
    otherwise
        tol = Inf;
        warning( 'CVX:SolverWarning', 'Gurobi returned an unknown status "%s". Please contact CVX Research Support.', res.status );
end
if isempty(x)
    x = NaN * ones(n,1); 
end
if need_y && isempty(y)
    y = NaN * ones(m,1); 
end
if need_z && isempty(z)
    z = NaN * ones(n,1); 
end
if tol == Inf
    status = 'Failed';
elseif tol > prec(2)
    status = [ 'Inaccurate/', status ];
end
if ~isempty( lbound )
    tol(2) = lbound;
end
iters = 0;

