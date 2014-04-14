function y = pos( x )

% POS    Positive part.
%    POS(X) = MAX(X,0). X must be real.
%
%    Disciplined convex programming information:
%        POS(X) is convex and nondecreasing in X. Thus when used in CVX
%        expressions, X must be convex (or affine).

persistent params
if isempty( params ),
    params.map = cvx_remap( ...
        { 'real' }, ...
        { 'p_nonconst' }, ...
        { 'n_nonconst' }, ...
        { 'r_affine', 'convex' } );
    params.constant = 1;
    params.funcs = { @pos_1, @pos_2, @pos_3, @pos_4 };
    params.name = 'pos';
end

try
    y = cvx_unary_op( params, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = pos_1( x )
y = max( x, 0 );

function y = pos_2( x )
y = x;

function y = pos_3( x )
y = zeros(size(x));

function y = pos_4( x ) %#ok
cvx_begin
    epigraph variable y(size(x)) nonnegative
    x <= y; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.