function y = log( x )

%   Disciplined convex programming information:
%       LOG(X) is concave and nondecreasing in X. When used in CVX
%       expressions, X must be concave.
%
%   Disciplined geometric programming information:
%       LOG(X) is typically not used in geometric programs. Technically it
%       possible to do so in certain advanced cases, because monomials and
%       posynomials are treated by CVX as log-affine and log-convex
%       constructs, respectively. However, such usage is undocumented and
%       will not be officially supported.

cvx_expert_check( 'log', x );

persistent P
if isempty( P ),
    P.map = cvx_remap( ...
        { 'nonpositive', 'n_nonconst' }, ...
        { 'positive', 'monomial' }, { 'posynomial' }, ...
        { 'concave' }, [0,2,3,4] );
    P.funcs = { [], @log_lkup, @log_posy, @log_gen };
end

try
    y = cvx_unary_op( P, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = log_lkup( x )
y = cvx( size( x ), cvx_getlog( cvx_basis( x ) ) );

function y = log_posy( x )
nx = numel( x );
x = cvx_basis( x );
[ rx, cx, vx ] = find( x );
nq = length( vx );
x = cvx( nq, sparse( rx, 1 : nq, vx ) );
cvx_begin gp
    variables w( nq )
    epigraph variable y( nx )
    exp( x ./ cvx_fastref( y, cx ) ) <= w; %#ok
    sparse( cx, 1 : nq, 1, nx, nq ) * log( w ) == 1; %#ok
cvx_end
y = log( y );

function y = log_gen( x ) %#ok
cvx_begin
    hypograph variable y(size(x))
    exp(y) <= x; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
