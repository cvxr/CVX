function y = geo_mean( x, dim, w )
error( nargchk( 1, 3, nargin ) );

%GEO_MEAN   Internal cvx version.

%
% Basic argument check
%

sx = size( x );
if nargin < 2 || isempty( dim ),
    dim = cvx_default_dimension( sx );
elseif ~cvx_check_dimension( dim ),
    error( 'Second argument must be a positive integer.' );
end

%
% Determine sizes, quick exit for empty arrays
%

sx = [ sx, ones( 1, dim - length( sx ) ) ];
nx = sx( dim );
sy = sx;
sy( dim ) = 1;
if any( sx == 0 ),
    y = ones( sy );
    return
end

%
% Third and fourth argument check
%

if nargin < 3 || isempty( w ),
    w = [];
elseif numel( w ) ~= length( w ) || ~isnumeric( w ) || ~isreal( w ) || any( w < 0 ) || any( w ~= floor( w ) ),
    error( 'Third argument must be a vector of nonnegative integers.' );
elseif length( w ) ~= nx,
    error( 'Third argument must be a vector of length %d', nx );
elseif ~any( w ),
    y = ones( sy );
    return
end

%
% Type check
%

persistent remap_1 remap_2 remap_3 remap_4
if isempty( remap_4 ),
    % Constant (postive or negative)
    remap_1 = cvx_remap( 'real' );
    remap_2 = cvx_remap( 'concave' );
    remap_3 = cvx_remap( 'log-convex' );
    remap_4 = cvx_remap( 'log-concave' );
end
vx = cvx_reshape( cvx_classify( x ), sx );
t1 = all( reshape( remap_1( vx ), sx ), dim );
t2 = all( reshape( remap_2( vx ), sx ), dim );
t3 = all( reshape( remap_3( vx ), sx ), dim ) | ...
     all( reshape( remap_4( vx ), sx ), dim );
% Valid combinations with zero or negative entries can be treated as constants
t1 = t1 | ( ( t2 | t3 ) & any( vx == 1 | vx == 9, dim ) );
ta = t1 + ( 2 * t2 + 3 * t3 ) .* ~t1;
nu = sort( ta(:) );
nu = nu([true;diff(nu)~=0]);
nk = length( nu );

%
% Permute and reshape, if needed
%

perm = [];
if nk > 1 || ( any( nu > 1 ) && nx > 1 ),
    if dim > 1 && any( sx( 1 : dim - 1 ) > 1 ),
        perm = [ dim, 1 : dim - 1, dim + 1 : length( sx ) ];
        x    = permute( x,  perm );
        sx   = sx( perm ); %#ok
        sy   = sy( perm );
        ta   = permute( ta, perm );
        dim  = 1;
    else
        perm = [];
    end
    nv = prod( sy );
    x  = reshape( x, nx, nv );
    ta = reshape( ta, 1, nv );
end

%
% Perform the computations
%

if nk > 1,
    y = cvx( [ 1, nv ], [] );
end
for k = 1 : nk,

    if nk == 1,
        xt = x;
        sz = sy; %#ok
    else
        tt = ta == nu( k );
        xt = cvx_subsref( x, ':', tt );
        nv = nnz( tt );
        sz = [ 1, nv ]; %#ok
    end

    switch nu( k ),
        case 0,
            error( 'Disciplined convex programming error:\n   Invalid computation: geo_mean( {%s} )', cvx_class( xt, true, true ) );
        case 1,
            yt = cvx( geo_mean( cvx_constant( xt ), dim, w ) );
        case 2,
            cvx_begin
                hypograph variable yt(sz);
                { cvx_accept_concave(xt), yt } == geo_mean_cone( size(xt), dim,  w, 'func' );
            cvx_end
        case 3,
            if nx == 1,
                yt = xt;
            elseif isempty( w ),
                yt = exp( sum( log( xt ), 1 ) * ( 1 / nx ) );
            else
                yt = exp( ( w / sum( w ) ) * log( xt ) );
            end
        otherwise,
            error( 'Shouldn''t be here.' );
    end

    if nk == 1,
        y = yt;
    else
        y = cvx_subsasgn( y, tt, yt );
    end

end

%
% Reverse the reshaping and permutation steps
%

y = reshape( y, sy );
if ~isempty( perm ),
    y = ipermute( y, perm );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
