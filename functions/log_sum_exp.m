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

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'real' ; 'convex' } );
    P.funcs = { @lse_1, @lse_2 };
    P.zero = -Inf;
    P.reduce = true;
    P.constant = 1;
    P.reverse = false;
    P.name = 'log_sum_exp';
    P.dimarg = 2;
end
cvx_expert_check( 'log_sum_exp', varargin{1} );
y = cvx_reduce_op( P, varargin{:} );

function x = lse_1( x )
xmid = 0.5 * ( max( x, [], 1 ) + min( x, [], 1 ) );
x = log( sum( exp( bsxfun( @minus, x, xmid ) ), 1 ) ) + xmid;

function y = lse_2( x ) %#ok
[nx,nv] = size(x); %#ok
cvx_begin
    variable w( nx, nv )
    epigraph variable y( 1, nv )
    exp( x - repmat( y, [ nx, 1 ] ) ) <= w; %#ok
    sum( w, 1 ) == 1; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

