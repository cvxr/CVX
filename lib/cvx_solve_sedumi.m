function [ x, y, status, tol, iters, z ] = cvx_solve_sedumi( At, b, c, nonls, quiet, prec )

n_in = 0;
n_out = 0;
max_sc = 0;
n = length( c );
m = length( b );
K = struct( 'l', 0, 'q', [], 's', [], 'scomplex', [], 'ycomplex', [] );
reord = struct( 'f', [], 'l', [], 'q', [], 'sr', [], 'sc', [], 'sv', [] );
for k = 1 : length( nonls ),
    temp = nonls( k ).indices;
    nn = size( temp, 1 );
    nv = size( temp, 2 );
    nnv = nn * nv;
    tt = nonls( k ).type;
    n_in = n_in + nnv;
    n_out = n_out + nnv;
    if nn == 1 || isequal( tt, 'nonnegative' ),
        K.l = K.l + nn * nv;
        reord.l = [ reord.l ; temp( : ) ];
    elseif isequal( tt, 'lorentz' ),
        K.q = [ K.q, nn * ones( 1, nv ) ];
        temp = temp( [ end, 1 : end - 1 ], : );
        reord.q = [ reord.q ; temp( : ) ];
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
        max_cc = max(cc);
        n_out = n_out + ( max_cc - nnv );
        cc = cc + max_sc;
        max_sc = max_sc + max_cc;
        reord.sr = [ reord.sr; rr( : ) ];
        reord.sc = [ reord.sc; cc( : ) ];
        reord.sv = [ reord.sv; vv( : ) ];
    end
end

n_f = n - n_in;
if n_f,
    n_out = n_out + n_f + 1;
    K.q(end+1) = n_f + 1;
    reord.f = ( 1 : n )';
    reord.f( [ reord.l ; reord.q ; reord.sr ] ) = [];
end
n_d = max( m - n - n_f + 1, isempty( At ) );
if n_d,
    n_out = n_out + n_d;
    K.l = K.l + n_d;
end
n_l = length(reord.l);
n_q = length(reord.q);
reord.sr = [ reord.l; reord.q;         reord.f;               reord.sr ];
reord.sc = [ 1:n_l,   n_l+n_d+(1:n_q), n_l+n_d+n_q+(2:n_f+1), reord.sc'+(n_out-max_sc) ];
reord.sv = [ ones(n_l+n_q+n_f,1);                             reord.sv ];
reord = sparse( reord.sr, reord.sc', reord.sv, n, n_out );

At = reord' * At;
c  = reord' * c;
pars.free   = 0;
pars.eps     = prec(1);
pars.bigeps  = prec(3);
if quiet,
    pars.fid = 0;
end
add_row = isempty( At );
if add_row,
    At = sparse( K.l, 1, 1, n_out, 1 );
    b = 1;
end

[ xx, yy, info ] = sedumi( At, b, c, K, pars );
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
else
    if isempty( status ),
        status = 'Solved';
    end
    if info.numerr == 1 && info.r0 > prec(2),
        status = [ 'Inaccurate/', status ];
    end
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
