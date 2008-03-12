function cvx_optval = logsumexp( x, dim )

%LOGSUMEXP    log(sum(exp(x))).
%   LOGSUMEXP(X) = LOG(SUM(EXP(X)). Unlike LOGSUMEXP_SDP, LOGSUMEXP is
%   computed using a successive approximation method, which provides
%   results exact to within the tolerance of the solver.
%
%   If X is a matrix, LOGSUMEXP_SDP(X) will perform its computations
%   along each column of X. If X is an N-D array, LOGSUMEXP_SDP(X)
%   will perform its computations along the first dimension of size
%   other than 1. LOGSUMEXP_SDP(X,DIM) will perform its computations
%   along dimension DIM.
%
%   Disciplined convex programming information:
%       LOGSUMEXP(X) is convex an nondecreasing in X; therefore, X
%       must be convex (or affine).

global cvx___
if ~cvx___.expert,
    error( sprintf( 'Disciplined convex programming error:\n    Exact logsumexp() distance is not yet supported.\n    Use LOGSUMEXP_SDP instead.' ) );
end

error( nargchk( 1, 2, nargin ) );
sx = size( x );
if nargin < 2 | isempty( dim ),
    dim = cvx_default_dimension( sx );
elseif ~cvx_check_dimension( dim, true ),
    error( 'Second argument must be a valid dimension.' );
end

%
% Quick exits
%

sx( end + 1 : dim ) = 1;
nx = sx( dim );
sy = sx;
sy( dim ) = 1;
if nx == 0,
    sx( dim ) = 1;
    cvx_optval = -Inf * ones( sy );
    return
elseif nx == 1,
    cvx_optval = x;
    return;
elseif any( sx == 0 ),
    cvx_optval = zeros( sy );
    return
end

%
% Determine the expression types
%

global cvx___
persistent remap
if isempty( remap ),
    remap_2 = cvx_remap( 'real' );
    remap_1 = cvx_remap( 'convex' ) & ~remap_2;
    remap = remap_1 + 2 * remap_2;
end
v = reshape( remap( cvx_classify( x ) ), sx );
v = min( v, [], dim );

%
% Process each type of expression one piece at a time
%

vu = unique( v );
nv = length( vu );
if nv > 1,
    y = cvx( sy, [] );
    if prod(sx(1:dim+1))>1 & prod(sx(dim+1:end))>1,
        perm = [ dim, 1:dim-1, dim+1:length(sx) ];
        x  = permute( x, perm );
        v  = permute( v, perm );
        y  = permute( y, perm );
        dim = 1;
    end
end
for k = 1 : nv,

    %
    % Select the category of expression to compute
    %

    vk = vu( k );
    if nv == 1,
        xt = x;
        sz = sy;
    else
        t = v == vk;
        xt = cvx_subsref( x, cvx_expand_dim( t, dim, nx ) );
        sz = size( xt );
        sx( dim ) = 1;
    end

    %
    % Perform the computations
    %

    switch vk,
        case 0,
            % Invalid
            error( sprintf( 'Disciplined convex programming error:\n    Illegal operation: logsumexp( {%s} ).', cvx_class( xt ) ) );
        case 1,
            % Affine, convex
            cvx_begin
                epigraph variable z( sz )
                sum( exp( cvx_accept_convex( x ) - cvx_expand_dim( z, dim, nx ) ), dim ) <= 1;
            cvx_end
        case 2,
            % Constant
            cvx_optval = logsumexp( cvx_constant( xt ) );
        otherwise,
            error( 'Shouldn''t be here.' );
    end

    %
    % Store the results
    %

    if nv == 1,
        y = cvx_optval;
    else
        y = cvx_subsasgn( y, t, cvx_optval );
    end

end

% Reshape again, just in case
y = reshape( y, sy );

% Copyright 2008 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
