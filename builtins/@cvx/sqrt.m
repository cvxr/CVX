function y = sqrt( x )

%   Discipined convex programming information for SQRT:
%      SQRT(X) is log-concave and nondecreasing in X. Therefore, when used
%      in DCPs, X must be concave (or affine).
%   
%   Disciplined geometric programming information for SQRT:
%      SQRT(X) is log-log-affine and nondecreasing in X. Therefore, when
%      used in DGPs, X may be log-affine, log-convex, or log-concave.

persistent remap funcs
if isempty( remap ),
    remap = cvx_remap( ...
        { 'constant' }, ...
        { 'l_valid' }, ...
        { 'r_affine', 'concave' } );
    funcs = { @sqrt_1, @sqrt_2, @sqrt_3 };
end

try
     y = unary_op( 'sqrt', funcs, remap, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = sqrt_1( x )
y = cvx( builtin( 'sqrt', cvx_constant( x ) ) );

function y = sqrt_2( x )
y = exp( 0.5 * log( x ) );

function y = sqrt_3( x )
sx = x.size_;
cvx_begin
    hypograph variable y( sx );
    square( x ) <= y; %#ok
    cvx_setnneg( y )
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
