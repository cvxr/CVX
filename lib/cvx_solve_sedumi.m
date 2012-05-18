function [ x, y, status, tol, iters, z ] = cvx_solve_sedumi( At, b, c, nonls, quiet, prec )

max_sc = 0;
n = length( c );
m = length( b );
K = struct( 'f', 0, 'l', 0, 'q', [], 's', [], 'scomplex', [], 'ycomplex', [] );

reord = struct( 'n', 0, 'r', [], 'c', [], 'v', [] );
reord = struct( 'f', reord, 'l', reord, 'a', reord, 'q', reord, 's', reord, 'h', reord );
reord.f.n = n;
for k = 1 : length( nonls ),
    temp = nonls( k ).indices;
    nn = size( temp, 1 );
    nv = size( temp, 2 );
    nnv = nn * nv;
    tt = nonls( k ).type;
    reord.f.n = reord.f.n - nnv;
    if nn == 1 || isequal( tt, 'nonnegative' ),
        reord.l.r = [ reord.l.r ; temp(:) ];
        reord.l.c = [ reord.l.c ; reord.l.n + ( 1 : nnv )' ];
        reord.l.v = [ reord.l.v ; ones(nnv,1) ];
        reord.l.n = reord.l.n + nnv;
    elseif isequal( tt, 'lorentz' ),
        if nn == 2,
            rr = [ temp ; temp ];
            cc = reshape( floor( 1 : 0.5 : 2 * nv + 0.5 ), 4, nv );
            vv = [1;1;-1;1]; vv = vv(:,ones(1,nv));
            reord.a.r = [ reord.a.r ; rr(:) ];
            reord.a.c = [ reord.a.c ; cc(:) + reord.a.n ];
            reord.a.v = [ reord.a.v ; vv(:) ];
            reord.a.n = reord.a.n + nnv;
        else
            temp = temp( [ end, 1 : end - 1 ], : );
            reord.q.r = [ reord.q.r ; temp(:) ];
            reord.q.c = [ reord.q.c ; reord.q.n + ( 1 : nnv )' ];
            reord.q.v = [ reord.q.v ; ones(nnv,1) ];
            reord.q.n = reord.q.n + nnv;
            K.q = [ K.q, nn * ones( 1, nv ) ];
        end
    else
        if isequal( tt, 'semidefinite' ),
            nn = 0.5 * ( sqrt( 8 * nn + 1 ) - 1 );
            str = cvx_create_structure( [ nn, nn, nv ], 'symmetric' );
        elseif isequal( tt, 'hermitian-semidefinite' ),
            K.scomplex = [ K.scomplex, length( K.s ) + ( 1 : nv ) ];
            nn = sqrt( nn );
            str = cvx_create_structure( [ nn, nn, nv ], 'hermitian' );
        else
            error( 'Unsupported nonlinearity: %s', tt );
        end
        K.s = [ K.s, nn * ones( 1, nv ) ];
        stri = cvx_invert_structure( str, 'compact' )';
        [ rr, cc, vv ] = find( stri );
        rr = temp( rr );
        reord.s.r = [ reord.s.r; rr( : ) ];
        reord.s.c = [ reord.s.c; cc( : ) + reord.s.n ];
        reord.s.v = [ reord.s.v; vv( : ) ];
        reord.s.n = reord.s.n + size( stri, 2 );
    end
end

if reord.f.n > 0,
    reord.f.r = [ 1 : n ]';
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
c  = reord' * c;
pars.free = K.f > 1 && nnz( K.q );
pars.eps     = prec(1);
pars.bigeps  = prec(3);
if quiet,
    pars.fid = 0;
end
add_row = isempty( At );
if add_row,
    K.f = K.f + 1;
    At = sparse( 1, 1, 1, n_out + 1, 1 );
    b = 1;
    c = [ 0 ; c ];
end
[ xx, yy, info ] = sedumi( At, b, c, K, pars );
if add_row,
    xx = xx(2:end);
    yy = zeros(0,1);
    At = zeros(n_out,0);
    b  = zeros(0,1);
    c  = c(2:end);
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
    z = - real( reord * ( At * yy ) );
    if add_row, y = zeros( 0, 1 ); end
elseif info.dinf ~= 0,
    status = 'Unbounded';
    y = NaN * ones( m, 1 );
    z = NaN * ones( n, 1 );
    x = real( reord * xx );
else
    x = real( reord * xx );
    y = yy;
    z = real( reord * ( c - At * yy ) );
    if add_row, y = zeros( 0, 1 ); end
end
if info.numerr == 2,
    status = 'Failed';
    if any( K.q == 2 ),
        warning( sprintf( ...
            [ 'This solver failure may possibly be due to a known bug in the SeDuMi solver.\n', ...
              'Try switching to SDPT3 by inserting "cvx_solver sdpt3" into your model.' ] ) );
    end
else
    if isempty( status ),
        status = 'Solved';
    end
    if info.numerr == 1 && info.r0 > prec(2),
        status = [ 'Inaccurate/', status ];
    end
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
