function shim = cvx_ecos( shim )

if ~isempty( shim.solve ),
    return
end
if isempty( shim.name ),
    fname = [ 'ecos.', mexext ];
    flen = length(fname);
    shim.name = 'ECOS';
    shim.config = struct( 'dualize', 1, 'nonnegative', 1, ...
        'lorentz', 1, 'exponential', 1 );
    fpaths = which( fname, '-all' );
    if ~iscell(fpaths),
        fpaths = { fpaths };
    end
    old_dir = pwd;
    oshim = shim;
    shim = [];
    for k = 1 : length(fpaths),
        fpath = fpaths{k};
        if ~exist( fpath, 'file' ) || any( strcmp( fpath, fpaths(1:k-1) ) ),
            continue
        end
        new_dir = fpath(1:end-flen-1);
        tshim = oshim;
        tshim.fullpath = fpath;
        tshim.version = 'unknown';
        if ~exist('OCTAVE_VERSION','builtin'),
            cd( new_dir );
            outp = evalc('ecos','[]'); %#ok
            [tok,remain] = strtok(outp,' ' ); %#ok
            tok = strtok(remain,' ' );
            if ~isempty(tok),
                tshim.version = tok;
            end
        end
        tshim.location = new_dir;
        if isempty( tshim.error ),
            tshim.solve = @solve;
            if k ~= 1,
                tshim.path = [ new_dir, pathsep ];
            end
        end
        shim = [ shim, tshim ]; %#ok
    end
    cd( old_dir );
    if isempty( shim ),
        shim = oshim;
        shim.error = 'https://github.com/ifa-ethz/ecos';
    end
else
    shim.solve = @solve;
end

function [ x, status, tol, iters, y, rawsol ] = solve( At, b, c, nonls, params )

USE_PRIMAL = 1;
rawsol = [];
n = length( c );
m = length( b );
K = struct( 'f', 0, 'l', 0, 'q', [], 'e', 0 );
reord = struct( 'n', 0, 'r', [], 'c', [], 'v', [] );
reord = struct( 'f', reord, 'l', reord, 'q', reord, 'e', reord );
reord.f.n = n;
for k = 1 : length( nonls ),
    temp = nonls( k ).indices;
    nn = size( temp, 1 );
    nv = size( temp, 2 );
    nnv = nn * nv;
    tt = nonls( k ).type;
    reord.f.n = reord.f.n - nnv;
    switch tt,
    case 'nonnegative',
        reord.l.r = [ reord.l.r ; temp(:) ];
        reord.l.c = [ reord.l.c ; reord.l.n + ( 1 : nnv )' ];
        reord.l.v = [ reord.l.v ; ones( nnv, 1 ) ];
        reord.l.n = reord.l.n + nnv;
    case 'lorentz',
        temp = temp( [ end, 1 : end - 1 ], : );
        reord.q.r = [ reord.q.r ; temp(:) ];
        reord.q.c = [ reord.q.c ; reord.q.n + ( 1 : nnv )' ];
        reord.q.v = [ reord.q.v ; ones(nnv,1) ];
        reord.q.n = reord.q.n + nnv;
        K.q = [ K.q, nn * ones( 1, nv ) ];
    case 'exponential',
        if USE_PRIMAL
            temp = temp([1,3,2],:);
        else
            temp = temp([3,1,2],:);
        end
        reord.e.r = [ reord.e.r ; temp(:) ];
        reord.e.c = [ reord.e.c ; reord.e.n + ( 1 : nnv )' ];
        if USE_PRIMAL,
            vv = ones(nnv,1);
        else
            vv = [-1;-1;exp(1)] * ones(1,nv);
        end
        reord.e.v = [ reord.e.v ; vv(:) ];
        reord.e.n = reord.e.n + nnv;
        K.e = K.e + nv;
    end
end
if reord.f.n > 0,
    reord.f.r = ( 1 : n )';
    reord.f.r( [ reord.l.r ; reord.q.r ; reord.e.r ] ) = [];
    reord.f.c = ( 1 : reord.f.n )';
    reord.f.v = ones(reord.f.n,1);
end
n_d = max( m - n - reord.f.n + 1, isempty( At ) );
if n_d,
    reord.l.n = reord.l.n + n_d;
end
K.f = reord.f.n;
K.l = reord.l.n;
n_out = reord.f.n;
reord.l.c = reord.l.c + n_out; n_out = n_out + reord.l.n;
reord.q.c = reord.q.c + n_out; n_out = n_out + reord.q.n;
reord.e.c = reord.e.c + n_out; n_out = n_out + reord.e.n;
reord = sparse( ...
    [ reord.f.r ; reord.l.r ; reord.q.r ; reord.e.r ], ...
    [ reord.f.c ; reord.l.c ; reord.q.c ; reord.e.c ], ...
    [ reord.f.v ; reord.l.v ; reord.q.v ; reord.e.v ], ...
    n, n_out );

At = reord' * At;
c  = reord' * c;

prec = params.precision;
opts.abstol = prec(1);
opts.reltol = prec(1);
opts.feastol = prec(1);
opts.abstol_inacc = prec(3);
opts.reltol_inacc = prec(3);
opts.feastol_inacc = prec(3);
opts.verbose = 1;
opts.maxit = 100;
if params.quiet,
    opts.verbose = 0;
end

if USE_PRIMAL,
    n_pos = n_out - K.f;
    ecos_c = full(c);
    ecos_G = sparse(1:n_pos,K.f+1:n_out,-1);
    ecos_h = zeros(n_pos,1);
    ecos_A = At';
    ecos_b = full(b);
else
    ecos_c = -full(b); % PHLI: Empirically this needed to be negated to match SeDuMi output
    ecos_G = At((K.f+1):end,:);
    ecos_h = full(c((K.f+1):end));
    ecos_A = At(1:K.f,:);
    ecos_b = full(c(1:K.f));
end
K.q = K.q(:);

varnames = {'x', 'y', 'info', 's', 'z'};
if USE_PRIMAL,
    [ xx, yy, info, xK, zz ] = cvx_run_solver( @ecos, ecos_c, ecos_G, ecos_h, K, ecos_A, ecos_b, opts, varnames{:}, 7, params ); %#ok
    yy = -yy;
else
    [ yy, xf, info, zz, xK ] = cvx_run_solver( @ecos, ecos_c, ecos_G, ecos_h, K, ecos_A, ecos_b, opts, varnames{:}, 7, params ); %#ok
    xx = [xf;xK];
end

if ~isfield( info, 'r0' ) && info.pinf,
    info.r0 = 0;
    info.iter = 0;
    info.numerr = 0;
end
tol = info.r0;
iters = info.iter;
xx = full( xx );
yy = full( yy );
status = '';
if info.pinf ~= 0,
    status = 'Infeasible';
    x = NaN * ones( n, 1 );
    y = yy;
elseif info.dinf ~= 0
    status = 'Unbounded';
    y = NaN * ones( m, 1 );
    x = real( reord * xx );
else
    x = real( reord * xx );
    y = yy;
end
if info.numerr == 2,
    status = 'Failed';
else
    if isempty( status ),
        status = 'Solved';
    end
    if info.numerr == 1 && info.r0 > prec(2),
        status = [ 'Inaccurate/', status ];
    end
end

% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.