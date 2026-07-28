function shim = cvx_glpk(shim)

% CVX_SOLVER_SHIM	GLPK interface for CVX.
%   This procedure returns a 'shim': a structure containing the necessary
%   information CVX needs to use this solver in its modeling framework.

global cvx___
if ~isempty(shim.solve)
    return
end

fs = cvx___.fs;
ps = cvx___.ps;

if isempty(shim.name)
    fname = 'glpk.m';
    flen = length(fname);
    int_path = [cvx___.where, fs];
    int_plen = length(int_path);
    shim.name = 'GLPK';
    shim.dualize = true;
    intname = [int_path, 'glpk', fs, fname];
    fpaths = which(fname, '-all');
    if exist(intname, 'file'), fpaths = [{intname}; fpaths]; end
    old_dir = pwd;
    oshim = shim;
    shim = [];
    for k = 1:length(fpaths)
        fpath = fpaths{k};
        if ~exist(fpath, 'file') || any(strcmp(fpath, fpaths(1:k-1)))
            continue
        end
        new_dir = fpath(1:end-flen-1);
        cd(new_dir);
        tshim = oshim;
        tshim.fullpath = fpath;
        tshim.version = 'unknown';
        is_internal = strncmp(new_dir, int_path, int_plen);
        if is_internal
            tshim.location = ['{cvx}', new_dir(int_plen:end)];
        else
            tshim.location = new_dir;
        end
        try
            fid = fopen(fname);
            otp = fread(fid, Inf, 'uint8=>char')';
            fclose(fid);
        catch errmsg
            tshim.error = sprintf('Unexpected error:\n%s\n', errmsg.message);
        end
        if isempty(tshim.error)
            otp = regexp(otp, 'Version: \d+\.\d+', 'match');
            if ~isempty(otp), tshim.version = otp{1}(10:end); end
            tshim.path = [new_dir, ps];
            tshim.check = @check;
            tshim.solve = @solve;
            tshim.eargs = {};
        end
        shim = [shim, tshim]; %#ok
    end
    cd(old_dir);
    if isempty(shim)
        shim = oshim;
        shim.error = 'Could not find a GLPK installation.';
    end
else
    shim.check = @check;
    shim.solve = @solve;
    shim.eargs = {};
end

function found_bad = check(nonls) %#ok
found_bad = false;

function [x, status, tol, iters, y, z] = solve(At, b, c, nonls, quiet, prec, settings)

n = length(c);
m  = length( b );
lb = -Inf(n, 1);
ub = +Inf(n, 1);
vtype = 'C';
vtype = vtype(ones(n, 1));
ctype = 'S';
ctype = ctype(ones(m, 1));
rr = zeros(0, 1);
cc = rr;
vv = rr;
zinv = rr;
is_ip = false;
for k = 1:length(nonls)
    temp = nonls(k).indices;
    nn = size(temp, 1);
    nv = size(temp, 2);
    tt = nonls(k).type;
    if strncmp(tt, 'i_', 2)
        is_ip = true;
        vartype(temp) = 'I';
        if strcmp(tt, 'i_binary')
            lb(temp) = 0;
            ub(temp) = 1;
        end
    elseif nn == 1 || isequal(tt, 'nonnegative')
        lb(temp) = 0;
    elseif isequal(tt, 'lorentz')
        if nn == 2
            rr2 = [temp; temp];
            cc2 = reshape(floor(1:0.5:2*nv+0.5), 4, nv);
            vv2 = [1; 1; -1; 1]; vv = vv(:, ones(1, nv));
            rr = [rr; rr(:)];
            cc = [cc; cc(:)];
            vv = [vv; vv(:)];
            zinv = [zinv; temp(:)];
        else
            error('GLPK does not support nonlinear constraints.');
        end
    else
        error('GLPK does not support nonlinear constraints.');
    end
end
if ~isempty(rr)
    znorm = [1:n]';
    znorm(zinv) = [];
    rr = [rr; znorm];
    cc = [cc; znorm];
    vv = [vv; ones(size(znorm))];
    reord = sparse(rr, cc, vv, n, n);
    At = reord'*At;
    c = reord'*c;
end
if quiet
    param.msglev = 0;
else
    param.msglev = 2;
end
param.scale = 128;
param.tolbnd = prec(1);
param.toldj = prec(1);
param.tolobj = prec(1);

[xx, fmin, errnum, extra] = cvx_run_solver(@glpk, c, At', b, lb, ub, ctype, vtype, 1, param, 'xx', 'fmin', 'errnum', 'extra', settings, 9);

tol = [];
iters = [];
x = full(xx);
y = full(extra.lambda);
z = full(extra.redcosts);
if ~isempty(rr)
    x = reord*x;
    z = reord*z;
    z(zinv) = z(zinv)*0.5;
end
status = 'Failed';
switch errnum
    case 0
        switch extra.status
            case 2
                if is_ip
                    status = 'Suboptimal';
                elseif errnumn == 0
                    status = 'Solved';
                else
                    status = 'Inaccurate/Solved';
                end
            case {3, 4}
                status = 'Infeasible';
            case 5
                status = 'Solved';
            case 6
                status = 'Unbounded';
        end
    case 10
        status = 'Infeasible';
    case 11
        status = 'Unbounded';
    case {5, 17, 15, 19}
        status = 'Failed';
    case {6, 7, 8, 9, 13, 14}
        switch extra.status
            case {2, 5}
                if is_ip
                    status = 'Suboptimal';
                else
                    status = 'Inaccurate/Solved';
                end
            case {3, 4}
                status = 'Inaccurate/Infeasible';
            case 6
                status = 'Inaccurate/Unbounded';
        end
end
if strcmp(status, 'Failed')
    tol = Inf;
elseif strncmp(status, 'Inaccurate/', 11)
    tol = prec(3);
else
    tol = prec(2);
end

% Copyright 2005-2016 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
