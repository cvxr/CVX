function [ outp, sdp_mode ] = cvx_pushcnstr( x, y, op, sdp_mode )

global cvx___
try
    pstr = cvx___.problems(end);
catch
    cvx_throw( 'No CVX model is present.' );
end

%
% Check problem
%

y_orig = y;
if isa( x, 'cvxcnst' ),
    x = rhs( x );
end

%
% Check for a dual reference
%

dx = cvx_getdual( x );
dy = cvx_getdual( y );
if isempty( dx ),
    dx = dy;
elseif ~isempty( dy ),
    cvx_throw( 'Two dual variable references found: %s %s', dx, dy );
end
if ~isempty( dx ),
    duals = pstr.duals;
    try
        dual = builtin( 'subsref', duals, dx );
    catch
        nm = cvx_subs2str( dx );
        cvx_throw( 'Dual variable "%s" has not been declared.', nm(2:end) );
    end
    if ~isempty( dual ),
        nm = cvx_subs2str( dx );
        cvx_throw( 'Dual variable "%s" is already in use.', nm(2:end) );
    end
end

%
% Check arguments
%

cx = isnumeric( x ) | isa( x, 'cvx' );
cy = isnumeric( y ) | isa( y, 'cvx' );
if ~cx,
    x = cvx_collapse( x, false, true );
    cx = isnumeric( x ) | isa( x, 'cvx' );
end
if ~cy,
    y = cvx_collapse( y, false, true );
    cy = isnumeric( y ) | isa( y, 'cvx' );
end
if ~cx || ~cy,
    if cx || cy || op(1) ~= '=',
        cvx_throw( 'Invalid CVX constraint: {%s} %s {%s}', class( x ), op, class( y ) );
    end
    sx = size( x );
    if ~isequal( sx, size( y ) ),
        cvx_throw( 'The left- and right-hand sides have incompatible sizes.' );
    else
        if ~isempty( dx ),
            pstr.duals = builtin( 'subsasgn', pstr.duals, dx, cell(sx) );
            cvx___.problems( end ) = pstr;
        end
        for k = 1 : prod( sx ),
            cvx_pushcnstr( x{k}, y{k}, op );
        end
        if nargout,
            outp = cvxcnst( y_orig );
        end
        return
    end
end

persistent param_eq param_ge param_le
if isempty( param_eq ),
    param_eq.map = cvx_remap( ...
            { { 'affine',  } }, ...
            { { 'l_affine' } } );
    param_eq.funcs = { @minus_nc, @eq_2 };
    param_eq.name = '=';
    param_eq.constant = [];
    param_ge.map = cvx_remap( ...
            { { 'constant' } }, ...
            { { 'l_concave' }, { 'l_convex' } }, ...
            { { 'concave' }, { 'convex' } }, ...
            [ 1, 2, 1 ] );
    param_ge.funcs = { @minus_nc, @ge_2 };
    param_ge.name = '>=';
    param_ge.constant = false;
    param_le.map = param_ge.map';
    param_le.funcs = { @le_1, @le_2 };
    param_le.name = '<=';
    param_le.constant = false;
end

if nargin < 5,
    sdp_mode = pstr.sdp;
end
switch op(1),
case '=',
    P = param_eq;
    sdp_mode = false;
case '<',
    P = param_le;
case '>',
    P = param_ge;
end
P.name = op;
P.sdp = sdp_mode;
[ z, sdp_mode ] = cvx_binary_op( P, x, y );
    
%
% Handle LMIs
%

zR = [];
sz = size( z );
zr = isreal( z );
if sdp_mode,
    z  = cvx_basis( z );
    dg = false;
    if zr,
        n  = sz(1);
        zq = any( z, 1 );
        nq = nnz( zq );
        if nq <= n,
            qn = bsxfun( @plus, (1:n+1:n*n)', 0:n*n:prod(sz)-1 );
            dg = nq == nnz( zq( qn ) );
            if dg,
                [ tx, dummy ] = find( cvx_basis( nonnegative( nq ) ) ); %#ok
                zz = sparse( tx, qn(zq(qn)), 1, tx(end), prod(sz) );
            end
        end
    end
    if ~dg,
        zt = z(:,permute(reshape(1:prod(sz),sz),[2,1,3:length(sz)]));
        if ~zr, zt = conj( zt ); end
        err = z - zt;
        if nnz(err),
            z   = 0.5 * ( z + zt );
            err = max( sum(abs(err),1) ./ max(realmin,sum(abs(z),1)) );
            if err > 8 * eps,
                if isreal(z), str = 'symmetric'; else str = 'Hermitian'; end
                cvx_throw( [ ...
                    'SDP constraints are expected to be %s (deviation: %g).\n', ...
                    '--- If this number is small (<1e-6), it may simply be due to roundoff error.\n', ...
                    '    This can be corrected by applying the SYM(X) function.\n', ...
                    '--- Otherwise, this is likely due to a modeling error. Did you declare the\n', ...
                    '    relevant matrix variables to be "symmetric" or "hermitian"?' ], str, full(err) );
            end
        end
        zz = cvx_basis( semidefinite( sz, ~zr ) );
    end
    z(size(zz,1),:) = 0;
    z = z - zz;
    iseq = true;
    if ~dg,
        [ zR, z ] = cvx_bcompress( z, 'full' );
    end
else
    iseq = op(1) == '=';
    z = cvx_basis( z );
end    

%
% Add the (in)equalities, create the dual variable if request
%

need_dual = ~isempty( dx );
zD = cvx_newcnstr( z, iseq, need_dual );
if need_dual
    if isempty( cvx_getdual( z ) ),
        if ~isempty( zR ), zD = zD * cvx_invert_structure( zR ); end
        pstr.duals = builtin( 'subsasgn', pstr.duals, dx, cvx( sz, zD ) );
        cvx___.problems( end ) = pstr;
    else
        warning( 'CVX:NoDualVariables', 'Dual variables are not availble for log-convex constraints.' )
    end
end

%
% Create the output object
%

if nargout,
    outp = cvxcnst( y_orig );
end

function z = eq_2( x, y )
z = cvx_setdual( minus_nc( log( x ), log( y ) ), 1 );

function z = ge_2( x, y )
z = cvx_setdual( minus_nc( log( x ), log( y ) ), 1 );

function z = le_1( x, y )
z = minus_nc( y, x );

function z = le_2( x, y )
z = cvx_setdual( minus_nc( log( y ), log( x ) ), 1 );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
