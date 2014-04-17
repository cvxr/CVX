function [ outp, sdp_mode ] = newcnstr( prob, x, y, op, sdp_mode )
    
%
% Check problem
%

if ~isa( prob, 'cvxprob' ),
    error( 'A cvx problem must be created first.' );
end
global cvx___
p = prob.index_;
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
    error( [ 'Two dual variable references found: "', dx, '","', dy, '"' ] );
end
if ~isempty( dx ),
    duals = cvx___.problems( p ).duals;
    try
        dual = builtin( 'subsref', duals, dx );
    catch
        nm = cvx_subs2str( dx );
        error( [ 'Dual variable "', nm(2:end), '" has not been declared.' ] );
    end
    if ~isempty( dual ),
        nm = cvx_subs2str( dx );
        error( [ 'Dual variable "', nm(2:end), '" already in use.' ] );
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
        error( 'Invalid CVX constraint: {%s} %s {%s}', class( x ), op, class( y ) );
    end
    sx = size( x );
    if ~isequal( sx, size( y ) ),
        error( 'The left- and right-hand sides have incompatible sizes.' );
    else
        if ~isempty( dx ),
            duals = cvx___.problems( p ).duals;
            duals = builtin( 'subsasgn', duals, dx, cell(sx) );
            cvx___.problems( p ).duals = duals;
        end
        for k = 1 : prod( sx ),
            newcnstr( prob, x{k}, y{k}, op );
        end
        if nargout,
            outp = cvxcnst( prob, y_orig );
        end
        return
    end
end

%
% Check readlevel
%

tx = cvx_readlevel( x );
ty = cvx_readlevel( y );
tx = any( tx( : ) > p );
ty = any( ty( : ) > p );
if tx || ty,
    error( 'Constraints may not involve internal, read-only variables.' );
end

persistent param_eq param_ge param_le
if isempty( param_eq ),
    param_eq.map = cvx_remap( ...
            { { 'affine',  } }, ...
            { { 'l_affine' } } );
    param_eq.funcs = { @minus_nc, @eq_2 };
    param_eq.name = '=';
    param_ge.map = cvx_remap( ...
            { { 'constant' } }, ...
            { { 'l_concave' }, { 'l_convex' } }, ...
            { { 'concave' }, { 'convex' } }, ...
            [ 1, 2, 1 ] );
    param_ge.funcs = { @minus_nc, @ge_2 };
    param_ge.name = '>=';
    param_le.map = param_ge.map';
    param_le.funcs = { @le_1, @le_2 };
    param_le.name = '<=';
end

try
    if nargin < 5,
        sdp_mode = cvx___.problems( p ).sdp;
    end
    switch op(1),
    case '=',
        params = param_eq;
        sdp_mode = false;
    case '<',
        params = param_le;
    case '>',
        params = param_ge;
    end
    params.name = op;
    params.sdp = sdp_mode;
    [ z, sdp_mode ] = cvx_binary_op( params, x, y );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc );
    else rethrow( exc ); end
end
    
%
% Handle LMIs
%

sz = size( z );
if sdp_mode,
    z  = cvx_basis( z );
    zr = isreal( z );
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
        zt = z(:,permute(reshape(1:prod(sz),sz),[2,1,sz(3:end)]));
        if ~zr, zt = conj( zt ); end
        err = z - zt;
        z   = 0.5 * ( z + zt );
        err = max( sum(abs(err),1) ./ max(realmin,sum(abs(z),1)) );
        if err > 8 * eps,
            if isreal(z), str = 'symmetric'; else str = 'Hermitian'; end
            error( 'CVX:ArgError', [ ...
                'SDP constraints are expected to be %s (deviation: %g).\n', ...
                '--- If this number is small (<1e-6), it may simply be due to roundoff error.\n', ...
                '    This can be corrected by applying the SYM(X) function.\n', ...
                '--- Otherwise, this is likely due to a modeling error. Did you declare the\n', ...
                '    relevant matrix variables to be "symmetric" or "hermitian"?' ], str, full(err) );
        end
        zz = cvx_basis( semidefinite( sz, ~zr ) );
    end
    z(size(zz,1),:) = 0;
    z = cvx( sz, z - zz );
    op = '==';
end

%
% Eliminate lexical redundancies
%

if op( 1 ) == '=',
    cmode = 'full';
else
    cmode = 'magnitude';
end
[ zR, zL ] = bcompress( z, cmode );

%
% Add the (in)equalities
%

touch( prob, zL, op(1) == '=' );
mO = length( cvx___.equalities );
mN = length( zL );
cvx___.equalities = vertcat( cvx___.equalities, zL );
cvx___.needslack( end + 1 : end + mN, : ) = op( 1 ) ~= '=';

%
% Create the dual
%

if ~isempty( dx )
    if isempty( cvx_getdual( z ) ),
        zI = cvx_invert_structure( zR )';
        zI = sparse( mO + 1 : mO + mN, 1 : mN, 1 ) * zI;
        zI = cvx( sz, zI );
        duals = builtin( 'subsasgn', duals, dx, zI );
        cvx___.problems( p ).duals = duals;
    else
        warning( 'CVX:NoDualVariables', 'Dual variables are not availble for log-convex constraints.' )
    end
end

%
% Create the output object
%

if nargout,
    outp = cvxcnst( prob, y_orig );
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
