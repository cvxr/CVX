function shim = cvx_sedumi( shim )

% CVX_SOLVER_SHIM	SeDuMi interface for CVX.
%   This procedure returns a 'shim': a structure containing the necessary
%   information CVX needs to use this solver in its modeling framework.

global cvx___
if ~isempty( shim.solve ),
    return
end
if isempty( shim.name ),
    fname = 'sedumi.m';
    fs = cvx___.fs;
    ps = cvx___.ps;
    int_path = [ cvx___.where, fs ];
    int_plen = length( int_path );
    shim.name = 'SeDuMi';
    shim.config = struct( 'dualize', 1, 'nonnegative', 1, ...
        'lorentz', -1, 'semidefinite', -1, ...
        'hermitian_semidefinite', -1 );
    flen = length(fname);
    fpaths = { [ int_path, 'sedumi', fs, fname ] };
    fpaths = [ fpaths ; which( fname, '-all' ) ];
    old_dir = pwd;
    oshim = shim;
    shim = [];
    for k = 1 : length(fpaths),
        fpath = fpaths{k};
        if ~exist( fpath, 'file' ) || any( strcmp( fpath, fpaths(1:k-1) ) ),
            continue
        end
        new_dir = fpath(1:end-flen-1);
        cd( new_dir );
        tshim = oshim;
        tshim.fullpath = fpath;
        tshim.version = 'unknown';
        is_internal = strncmp( new_dir, int_path, int_plen );
        if is_internal,
            tshim.location = [ '{cvx}', new_dir(int_plen:end) ];
        else
            tshim.location = new_dir;
        end
        try
            fid = fopen(fname);
            otp = fread(fid,Inf,'uint8=>char')';
            fclose(fid);
        catch errmsg
            tshim.error = sprintf( 'Unexpected error reading %s:\n%s\n', fname, errmsg.message );
        end
        if isempty( tshim.error ),
            otp = regexp( otp, 'SeDuMi \d\S+', 'match' );
            if ~isempty(otp), tshim.version = otp{end}(8:end); end
            vnum = str2double( tshim.version );
            tshim.solve = @solve;
            tshim.config.convertSDP = vnum >= 1.3 && vnum < 1.32;
            if k ~= 2,
                tshim.path = [ new_dir, ps ];
                if ~isempty(cvx___.msub) && exist([new_dir,fs,cvx___.msub],'dir'),
                   tshim.path = [ new_dir, fs, cvx___.msub, ps, tshim.path ];
                end
            end
        end
        shim = [ shim, tshim ]; %#ok
    end
    cd( old_dir );
    if isempty( shim ),
        shim = oshim;
    end
else
    shim.solve = @solve;
end

function [ x, status, tol, iters, y ] = solve( At, b, c, nonls, params )

n = length( c );
m = length( b );
K = struct( 'f', 0, 'l', 0, 'q', [], 'r', [], 's', [], 'scomplex', [], 'ycomplex', [] );
reord = struct( 'n', 0, 'r', [], 'c', [], 'v', [] );
reord = struct( 'f', reord, 'l', reord, 'a', reord, 'q', reord, 's', reord );
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
        case 'semidefinite',
            nn = 0.5 * ( sqrt( 8 * nn + 1 ) - 1 );
            str = cvx_create_structure( [ nn, nn, nv ], 'symmetric' );
            K.s = [ K.s, nn * ones( 1, nv ) ];
            [ rr, cc, vv ] = find( cvx_invert_structure( str, true ) );
            rr = temp( rr );
            reord.s.r = [ reord.s.r; rr( : ) ];
            reord.s.c = [ reord.s.c; cc( : ) + reord.s.n ];
            reord.s.v = [ reord.s.v; vv( : ) ];
            reord.s.n = reord.s.n + nn * nn * nv;
        case 'hermitian_semidefinite',
            % SeDuMi's complex SDP support was restored in v1.33.
            K.scomplex = [ K.scomplex, length( K.s ) + ( 1 : nv ) ];
            nn = sqrt( nn );
            str = cvx_create_structure( [ nn, nn, nv ], 'hermitian' );
            K.s = [ K.s, nn * ones( 1, nv ) ];
            [ rr, cc, vv ] = find( cvx_invert_structure( str, true ) );
            rr = temp( rr );
            reord.s.r = [ reord.s.r; rr( : ) ];
            reord.s.c = [ reord.s.c; cc( : ) + reord.s.n ];
            reord.s.v = [ reord.s.v; vv( : ) ];
            reord.s.n = reord.s.n + nn * nn;
        otherwise,
            cvx_throw( 'Unexpected nonlinearity: %s', tt );
    end
end
if reord.f.n > 0,
    reord.f.r = ( 1 : n )';
    reord.f.r( [ reord.l.r ; reord.a.r ; reord.q.r ; reord.s.r ] ) = [];
    reord.f.c = ( 1 : reord.f.n )';
    reord.f.v = ones(reord.f.n,1);
end
n_d = max( m - n - reord.f.n + 1, isempty( At ) );
if n_d,
    reord.l.n = reord.l.n + n_d;
end
K.f = reord.f.n;
K.l = reord.l.n + reord.a.n;
n_out = reord.f.n;
reord.l.c = reord.l.c + n_out; n_out = n_out + reord.l.n;
reord.a.c = reord.a.c + n_out; n_out = n_out + reord.a.n;
reord.q.c = reord.q.c + n_out; n_out = n_out + reord.q.n;
reord.s.c = reord.s.c + n_out; n_out = n_out + reord.s.n;
reord = sparse( ...
    [ reord.f.r ; reord.l.r ; reord.a.r ; reord.q.r ; reord.s.r ], ...
    [ reord.f.c ; reord.l.c ; reord.a.c ; reord.q.c ; reord.s.c ], ...
    [ reord.f.v ; reord.l.v ; reord.a.v ; reord.q.v ; reord.s.v ], ...
    n, n_out );

At = reord' * At;
 c = reord' * c;
pars.free = K.f > 1 && nnz( K.q );
prec = params.precision;
pars.eps     = prec(1);
pars.bigeps  = prec(3);
if params.quiet,
    pars.fid = 0;
end
add_row = isempty( At );
if add_row,
    K.f = K.f + 1;
    At = sparse( 1, 1, 1, n_out + 1, 1 );
    b = 1;
    c = [ 0 ; c ];
end
[ xx, yy, info ] = cvx_run_solver( @sedumi, At, b, c, K, pars, 'xx', 'yy', 'info', 5, params );
if add_row,
    xx = xx(2:end);
    yy = zeros(0,1);
end
if ~isfield( info, 'r0' ) && info.pinf,
    info.r0 = 0;
    info.iter = 0;
    info.numerr = 0;
end
tol    = info.r0;
iters  = info.iter;
xx     = full( xx );
yy     = full( yy );
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

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
