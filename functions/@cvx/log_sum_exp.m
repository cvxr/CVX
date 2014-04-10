function y = log_sum_exp( varargin ) % x, dim

%LOG_SUM_EXP   CVX internal version.

cvx_expert_check( 'log_sum_exp', x );

persistent params
if isempty( params ),
    params.map = cvx_remap( { 'real' ; 'convex' } );
    params.funcs = { @lse_1, @lse_2 };
    params.zero = -Inf;
    params.reduce = true;
    params.reverse = false;
    params.dimarg = 2;
    params.name = 'log_sum_exp';
end

try
    y = reduce_op( params, varargin{:} );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end

function x = lse_1( x )
if x.size_(1) > 1,
    x = cvx_constant( x );
    xmid = 0.5 * ( max( x, [], 1 ) + min( x, [], 1 ) );
    x = cvx( log( sum( exp( bsxfun( @minus, x, xmid ) ) ) ) + xmid );
end

function y = lse_2( x )
if x.size_(1) > 1,
    w = []; y = [];
    [nx,nv] = size(x);
    cvx_begin
        variable w( nx, nv )
        epigraph variable y( 1, nv )
        { cvx_accept_convex( x ) - repmat( y, [nx,1] ), 1, w } == exponential( [nx,nv] ); %#ok
        sum( w, 1 ) == 1; %#ok
    cvx_end
else
    y = x;
end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
