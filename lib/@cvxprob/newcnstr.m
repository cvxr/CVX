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

if isa( x, 'cvxcnst' ), 
    x = rhs( x ); 
end
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

%
% Check sizes
%

dvalid = true;
if nargin < 5,
    sdp_mode = op(1) ~= '=' && cvx___.problems( p ).sdp;
end
[ x, y, sz, xs, ys, sdp_mode ] = cvx_broadcast( x, y, sdp_mode );
    
%
% Handle the SDP case
%

if sdp_mode,
    
    if op(1) == '>',
        z = plus( x, y, '-' );
    else
        z = plus( y, x, '-' );
    end
    zq = any( cvx_basis( z ), 1 );
    qn = bsxfun( @plus, (1:sz(1)+1:sz(1)*sz(2))', 0:sz(1)*sz(2):prod(sz)-1 );
    if nnz( zq ) == nnz( zq( qn ) ),
        [ tx, dummy ] = find( cvx_basis( nonnegative( numel(qn) ) ) ); %#ok
        zz = cvx( sz, sparse( tx, qn, 1, tx(end), prod(sz) ) );
    else
        zz = semidefinite( sz, ~isreal( x ) || ~isreal( y ) );
    end
    z = plus( z, zz, '-' );
    op = '==';

else
    
    persistent map_eq map_le map_ge map_ne map_N %#ok
    if isempty( map_ge ),
        map_ne = cvx_remap( ...
            { { } } );
        map_eq = cvx_remap( ...
            { { 'affine' } }, { { 'l_affine' } } );
        map_ge = cvx_remap( ...
            { { 'constant' } }, ...
            { { 'positive' }, { 'l_convex' } }, ...
            { { 'l_concave' }, { 'l_convex' } }, ...
            { { 'concave' }, { 'convex' } }, ...
            [ 1, 3, 2, 1 ] );
        map_le = map_ge';
        map_N  = size( map_ne, 1 );
    end
    switch op(1),
        case '<',
            remap = map_le;
        case '>',
            remap = map_ge;
        case '~',
            remap = map_ne;
        otherwise,
            remap = map_eq;
    end
    vx = cvx_classify( x );
    vy = cvx_classify( y );
    vm = vx + map_N * ( vy - 1 );
    vr = remap( vm );
    tt = vr == 0;
    if any( tt(:) ),
        cvx_dcp_error( x, y, tt, op );
    end
    % If the user wants to use a dual variable, we need to be more
    % conservative in judging which constraints we need to take the
    % logarithm of. Othewise, we can do them all
    if isempty( dx ),
        tt = vr >= 2;
    else
        tt = vr == 2;
    end
    if any( tt(:) ),
        dvalid = false;
        if ~isempty( dx ),
            tt = vr >= 2;
        end
        if all( tt(:) ),
            x = log( x );
            y = log( y );
        else
            if xs, x = x * ones(sz); end
            if ys, y = y * ones(sz); end
            x( tt ) = log( x( tt ) );
            y( tt ) = log( y( tt ) );
        end
    end
    if op(1) == '<',
        z = plus( y, x, '-' );
    else
        z = plus( x, y, '-' );
    end
    
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
if sdp_mode,
    if isreal( zR ),
        nnq = 0.5 * sz( 1 ) * ( sz( 1 ) + 1 );
    else
        nnq = sz( 1 ) * sz( 1 );
    end
    if size( zR, 1 ) > nnq * prod( sz( 3 : end ) ),
        warning( 'CVX:UnsymmetricLMI', [ ...
            'This linear matrix inequality appears to be unsymmetric. This is\n', ...
            'very likely an error that will produce unexpected results. Please check\n', ...
            'the LMI; and, if necessary, re-enter the model.' ], 1 );  
    end
end

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

if ~isempty( dx ),
    if dvalid,
        zI = cvx_invert_structure( zR )';
        zI = sparse( mO + 1 : mO + mN, 1 : mN, 1 ) * zI;
        zI = cvx( sz, zI );
        duals = builtin( 'subsasgn', duals, dx, zI );
        cvx___.problems( p ).duals = duals;
    else
        warning( 'CVX:NoDualVariables', 'Dual variables are not availble for this constraint.' )
    end
end

%
% Create the output object
%

if nargout,
    outp = cvxcnst( prob, y_orig );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
