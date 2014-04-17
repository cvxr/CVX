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
    P.map = cvx_remap( { 'positive' }, ...
        { 'l_concave_', 'l_affine', 'monomial' }, { 'l_convex' }, ...
        { 'r_affine', 'concave' } );
    P.funcs = { @log_cnst, @log_mono, @log_posy, @log_cncv };
end

try
    y = cvx_unary_op( P, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = log_cnst( x )
y = builtin( 'log', x );

function y = log_mono( x )
% Monomial
global cvx___
nb = prod( x.size_ );
[ rx, cx, vx ] = find( x.basis_ );
logs = cvx___.logarithm( rx, 1 );
tt = vx ~= 1; nt = sum( tt );
bx = sparse( [ ones( nt, 1 ) ; logs ], [ cx( tt ) ; cx ], [ log( vx( tt ) ) ; ones( nb, 1 ) ], ...
    full( max( logs ) ), size( x.basis_, 2 ) );
y = cvx( x.size_, bx );

function y = log_posy( x )
% Posynomial
global cvx___
x = x.basis_;
sx = x.size_;
rc = full( sum( x ~= 0, 1 ) );
ru = sort( rc(:) );
ru = ru([true;diff(ru)~=0]);
nu = length( ru );
if nu ~= 1,
    y = cvx( sx, [] );
end
for kk = 1 : nu,
    rk = ru( kk );
    if nu == 1,
        xt = x;
    else
        tt  = rc == rk;
        xt = x( :, tt );
    end
    [ rx, cx, vx ] = find( xt ); %#ok
    rx = rx( : ); vx = vx( : );
    nq = length( vx );
    vx = log( vx );
    tz = rx ~= 1;
    rx = cvx___.logarithm( rx( tz ), 1 );
    vx = vx + cvx( nq, sparse( rx, find( tz ), 1, full( max( rx ) ), nq ) );
    vx = reshape( vx, rk, nq / rk );
    vx = log_sum_exp( vx );
    if nu == 1,
        y = reshape( vx, sx );
    else
        y = cvx_subsasgn( y, tt, vx );
    end
end

function y = log_cncv( x ) %#ok
cvx_begin
    hypograph variable y( size(x) )
    exp( y ) <= x; %#ok
cvx_end


% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
