function y = exp( x )

%   Disciplined convex programming information:
%       EXP(X) is convex and nondecreasing in X. When used in CVX
%       expressions, X must be real. Typically, X must also be affine
%       or convex; X can also be concave, but this produces a log-concave
%       result with very limited usefulness.
%
%   Disciplined geometric programming information:
%       EXP(X) is typically not used in geometric programs. However,
%       EXP(X), where X is a monomial or posynomial, can be included in 
%       geometric programs wherever a posynomial would be appropriate.

cvx_expert_check( 'exp', x );

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'real' }, { 'convex', 'concave' } );
    P.funcs = { @exp_cnst, @exp_cxcv };
end

try
    y = cvx_unary_op( P, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = exp_cnst( x )
y = builtin( 'exp', x );

function y = exp_cxcv( x )
global cvx___
persistent expv ccv
if isempty( expv )
    expv = int8([3,3,3,4,15,15,15,16,16,16,17,17,17,19,19,17,17,17,19]);
end
x = sparsify( x, 'exponential' );
[ rx, cx, vx ] = find( x.basis_ );
tt = rx == 1;  rx( tt ) = [];
cc = cx( tt ); cx( tt ) = [];
vc = vx( tt ); % vx( tt ) = [];
exps = cvx___.exponential( rx, 1 );
tt = exps == 0;
if any( tt ),
    n1 = unique( rx( tt ) );
    n2 = newvar( cvx___.problems( end ).self, '', length( n1 ) );
    [ n2, dummy ] = find( n2.basis_ ); %#ok
    cvx___.exponential( n1, 1 ) = n2( : );
    cvx___.logarithm( n2, 1 ) = n1( : );
    cvx___.classes( n2 ) = expv( cvx___.classes( n1 ) );
    exps = cvx___.exponential( rx, 1 );
    cvx___.exp_used = true;
end
nb = size( x.basis_, 2 );
bx = sparse( exps, cx, 1, full( max( exps ) ), nb );
if ~isempty( cc ),
    bx = bx * diag(exp(sparse(cc,1,vc,nb,1)));
end
y = cvx( x.size_, bx );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
