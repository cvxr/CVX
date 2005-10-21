function [ value, x, y, status ] = cvx_solve_sedumi( A, b, c, d, nonls, quiet )
if nargin < 6, quiet = false; end

[ m, n ] = size( A );
x = zeros( n, 1 );
y = zeros( m, 1 );

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
K.l = K.l + K.a;
pars.eps = 1e-8;
pars.bigeps = 1e-3;
pars.cg.refine = 3;
if quiet,
    pars.fid = 0;
end
pars.free = K.f + K.l < n;

if isempty( A ),
    x = zeros( size( A, 2 ), 1 );
    y = zeros( 0, 1 );
    status = 'Solved';
    value = c' * x + d;
else,
    [ x, y, info ] = sedumi( A, b, c, K, pars );
    if info.pinf ~= 0,
        status = 'Infeasible';
        value = +Inf;
    elseif info.dinf ~= 0,
        status = 'Unbounded';
        value = -Inf;
    else,
        value = c' * x + d;
        if info.numerr == 0,
            status = 'Solved';
        else,
            if info.numerr == 1,
                status = 'Inaccurate';
            else
                status = 'Failed';
            end
            if info.feasratio > 0.99,
                status = [ status, ', likely close to a solution' ];
            elseif info.feasratio < -0.99,
                status = [ status, ', likely infasible' ];
                value = +Inf;
            end
        end
    end
end

x = real( reord * x );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
