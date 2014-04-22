function y = log_sum_exp( varargin )

%LOG_SUM_EXP    log(sum(exp(x))).
%   LOG_SUM_EXP(X) = LOG(SUM(EXP(X)).
%
%   When used in a CVX model, LOG_SUM_EXP(X) causes CVX's successive
%   approximation method to be invoked, producing results exact to within
%   the tolerance of the solver. This is in contrast to LOGSUMEXP_SDP,
%   which uses a single SDP-representable global approximation.
%
%   If X is a matrix, LOGSUMEXP_SDP(X) will perform its computations
%   along each column of X. If X is an N-D array, LOGSUMEXP_SDP(X)
%   will perform its computations along the first dimension of size
%   other than 1. LOGSUMEXP_SDP(X,DIM) will perform its computations
%   along dimension DIM.
%
%   Disciplined convex programming information:
%       LOG_SUM_EXP(X) is convex and nondecreasing in X; therefore, X
%       must be convex (or affine).

cvx_expert_check( 'log_sum_exp', x );

persistent params
if isempty( params ),
    params.map = cvx_remap( { 'real' ; 'convex' } );
    params.funcs = { @lse_1, @lse_2 };
    params.zero = -Inf;
    params.reduce = true;
    params.constant = 1;
    params.reverse = false;
    params.name = 'log_sum_exp';
    params.dimarg = 2;
end

try
    y = cvx_reduce_op( params, varargin{:} );
scatch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function x = lse_1( x )
xmid = 0.5 * ( max( x, [], 1 ) + min( x, [], 1 ) );
x = log( sum( exp( bsxfun( @minus, x, xmid ) ), 1 ) ) + xmid;

function y = lse_2( x ) %#ok
[nx,nv] = size(x);
cvx_begin
    variable w( nx, nv )
    epigraph variable y( 1, nv )
    { linearize( x ) - repmat( y, [nx,1] ), 1, w } == exponential( [nx,nv] ); %#ok
    sum( w, 1 ) == 1; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

