function [ value, x, y, status ] = cvx_solve_sedumi( A, b, c, d, sgn, nonls, quiet )
if nargin < 6, quiet = false; end

[ m, n ] = size( A );
K.f = 0; K.l = 0; K.a = 0; K.q = []; K.s = [];
reord.f = []; reord.l = []; reord.q = [];
reord.ar = []; reord.ac = []; reord.av = [];
reord.sr = []; reord.sc = []; reord.sv = [];
K.s = []; K.scomplex = []; K.ycomplex = [];
for k = 1 : length( nonls ),
    temp = nonls( k ).indices - 1;
    nn = size( temp, 1 );
    nv = size( temp, 2 );
    tt = nonls( k ).type;
    switch tt,
        case 'nonnegative',
            K.l = K.l + nn * nv;
            reord.l = [ reord.l ; temp( : ) ];
        case 'lorentz',
            if nn > 2,
                K.q = [ K.q, nn * ones( 1, nv ) ];
                temp = temp( [ end, 1 : end - 1 ], : );
                reord.q = [ reord.q ; temp( : ) ];
            else,
                K.a  = K.a + nn * nv;
                stri = cvx_replicate_structure( [ 1, -1 ; 1, 1 ], nv );
                [ rr, cc, vv ] = find( stri );
                if ~isempty( reord.ac ),
                    cc = cc + reord.ac( end );
                end
                rr = temp( rr );
                reord.ar = [ reord.ar; rr( : ) ];
                reord.ac = [ reord.ac; cc( : ) ];
                reord.av = [ reord.av; vv( : ) ];
            end
        case { 'semidefinite', 'hermitian-semidefinite' },
            if tt( 1 ) == 'h',
                K.scomplex = [ K.scomplex, length( K.s ) + ( 1 : nv ) ];
                nn = sqrt( nn );
                str = cvx_create_structure( [ nn, nn, nv ], 'hermitian' );
            else,
                nn = 0.5 * ( sqrt( 8 * nn + 1 ) - 1 );
                str = cvx_create_structure( [ nn, nn, nv ], 'symmetric' );
            end
            K.s = [ K.s, nn * ones( 1, nv ) ];
            stri = cvx_invert_structure( str, 'compact' );
            [ rr, cc, vv ] = find( stri );
            if ~isempty( reord.sc ),
                cc = cc + reord.sc( end );
            end
            rr = temp( rr );
            reord.sr = [ reord.sr; rr( : ) ];
            reord.sc = [ reord.sc; cc( : ) ];
            reord.sv = [ reord.sv; vv( : ) ];
        otherwise,
            error( sprintf( 'Nonlinearity "%s" not supported for SeDuMi', tt ) );
    end
end

reord.f = 1 : n;
reord.f( [ reord.l ; reord.ar ; reord.q ; reord.sr ] ) = [];
K.f = length( reord.f );
g1 = K.f + K.l;
g3 = sum( K.q );
g2 = K.f + K.l + K.a + g3;
reord.sr = [ reord.f(:); reord.l; reord.ar;    reord.q;       reord.sr    ];
reord.sc = [ [1:g1]';             reord.ac+g1; [g2-g3+1:g2]'; reord.sc+g2 ];
reord.sv = [ ones(g1,1);          reord.av;    ones(g3,1);    reord.sv    ];
reord = sparse( reord.sr, reord.sc, reord.sv );

A = A * reord;
c = reord' * c;
nn = size( A, 2 );
K.l = K.l + K.a;
pars.eps = 1e-8;
pars.bigeps = 1e-3;
pars.cg.refine = 3;
if quiet,
    pars.fid = 0;
end
pars.free = 0; % ~isempty( K.q ) | ~isempty( K.s );
pars.sdp = 0;

if ~quiet,
    disp( ' ' );
    disp( sprintf( 'Calling SeDuMi: %d variables (%d free), %d equality constraints', n, K.f, m ) );
    disp( '------------------------------------------------------------------------' );
end

if all( b == 0 ) & all( c == 0 ),
    if ~quiet,
        disp( 'Degenerate problem encountered; solution determined analytically.' );
    end
    status = 'Solved';
    x = zeros( nn, 1 );
    y = zeros( m, 1 );
    value = d;
elseif any( ~any( A, 2 ) ),
    if ~quiet,
        disp( 'Degenerate infeasible problem encountered; solution determined analytically.' );
    end
    status = 'Infeasible';
    value = -sgn * Inf;
    x = NaN * ones( nn, 1 );
    y = sgn * b .* ~any( A, 2 );
    y = y / abs( b' * y );
elseif isempty( b ),
    if ~quiet,
        disp( 'Degenerate problem encountered; an unbounded direction search will be' );
        disp( 'performed. The SeDuMi status messages will not coincide with cvx_status.' );
        disp( '------------------------------------------------------------------------' );
    end
    [ x, y, info ] = sedumi( sgn * c, -1, [], K, pars );
    x = full( x );
    y = zeros( 0, 1 );
    if info.numerr == 2 | info.dinf ~= 0,
        status = 'Failed';
        x = NaN * ones( nn, 1 );
        y = NaN * ones( m, 1 );
        value = NaN;
    elseif info.pinf ~= 0,
        status = 'Solved';
        x = zeros( nn, 1 );
        value = d;
    else,
        status = 'Unbounded';
        value = -Inf * sgn;
    end
    if info.numerr == 1 & pars.eps > 0,
        status = [ 'Inaccurate/', status ];
    end
else,
    degen = false;
    if K.f == size(A,2),
        degen = true;
        A( :, end + 1 ) = 0;
        K.l = 1;
    end
    [ x, y, info ] = sedumi( A, b, sgn * c, K, pars );
    x = full( x );
    y = full( y );
    if degen,
        x(end) = [];
    end
    if info.numerr == 2,
        status = 'Failed';
        x = NaN * ones( nn, 1 );
        y = NaN * ones( m, 1 );
        value = NaN;
    elseif info.pinf ~= 0,
        status = 'Infeasible';
        x = NaN * ones( nn, 1 );
        value = +Inf * sgn;
    elseif info.dinf ~= 0,
        status = 'Unbounded';
        y = NaN * ones( m, 1 );
        value = -Inf * sgn;
    elseif info.feasratio < 0,
        status = 'Infeasible';
        x = NaN * ones( nn, 1 );
        y = NaN * ones( m, 1 );
        value = +Inf * sgn;
    else,
        status = 'Solved';
        value = c' * x + d;
    end
    if info.numerr == 1 & pars.eps > 0,
        status = [ 'Inaccurate/', status ];
    end
end

if ~quiet,
    disp( '------------------------------------------------------------------------' );
    disp( sprintf( 'Optimal value (cvx_optval): %+g\nStatus (cvx_status): %s', full( value ), status ) );
    disp( ' ' );
end

x = real( reord * x );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
