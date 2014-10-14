function y = huber_pos( x, M, t )

%HUBER_POS   Monotonic huber penalty function.
%   For a vector X, HUBER_POS(X) computes the monotonic Huber-style function
% 
%                      0     if 0>=X
%       HUBER_POS(X) = X^2   if 0<=X<=1
%                      2*X-1 if    X>=1
%
%   HUBER_POS(X,M) is the monotonic Huber-style penalty function of
%   halfwidth M, M.^2.*HUBER_POS(X./M). M must be real and positive.
%
%   HUBER_POS(X,M,T) computes the monotonic Huber-style penalty function 
%   with halfwidth M and concomitant scale T:
%
%       HUBER_POS(X,M,T) = T.*HUBER_POS(X./T,M) if T > 0
%                          +Inf                 if T <= 0
%
%   See the help file for HUBER for information about this usage.
%
%   For matrices and N-D arrays, the penalty function is applied to each
%   element of X independently. M and T must be compatible with X in the same
%   sense as .*: one must be a scalar, or they must have identical size.
%
%   Disciplined convex programming information:
%       HUBER_POS is jointly convex in X and T. It is nondecreasing in X and
%       nonincreasing in T. Therefore, when used in CVX specifications, X
%       must be convex and T must be concave (or affine). Both must be real.

persistent P
if isempty( P ),
    P.map = cvx_remap( ...
        { { 'any' }, { 'nonpositive' } }, ...
        { { 'real' } }, ...
        { { 'convex' }, { 'concave' } }, [0,1,2] );
    P.funcs = { @huber_pos_c, @huber_pos_nc };
    P.constant = 1;
    P.name = 'huber_pos';
end
if nargin < 2,
    M = 1;
elseif ~( isnumeric(M) && numel(M)==1 && isreal(M) && M>0 ),
    cvx_throw( 'Second argument must be a positive scalar.' );
end
if nargin < 3,
    t = 1;
end
y = cvx_binary_op( P, x, t, M );

function z = huber_pos_c( x, t, M )
y = max( x, 0 );
z = min( y, M );
z = t .* z .* ( 2 * y - z );

function cvx_optval = huber_pos_nc( x, t, M ) %#ok
sz = max(numel(x),numel(t)); %#ok
cvx_begin
    variable w(sz) nonnegative
    variable v(sz) nonnegative
    minimize( quad_over_lin( w, t, 0 ) + ( 2 * M ) * v )
    x <= w + v; %#ok
    w <= M * t; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
