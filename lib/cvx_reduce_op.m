function y = cvx_reduce_op( p, varargin )

% Get dimension and main argument
if isempty( p.dimarg ),
    x = varargin{1};
    dim = varargin{2};
    zx = size( x );
    varargin([1,2]) = [];
else
    [ zx, x, dim, varargin ] = cvx_get_dimension( p.dimarg, varargin, 1 );
end
ox = isa( x, 'cvx' );
nd = length( zx );

% Quick exit for empty results
nx = zx( dim );
zy = zx;
if p.reduce && ( nx || ~isempty( p.zero ) ),
    zy( dim ) = 1; 
end
if ~all( zy ),
    y = cvx( repmat( p.zero ), zy );
    return
end
nv = prod( zx ) / nx;

% Permute if needed, and reshape to canonical size
ndxs = [];
if nx > 1 && nv > 1,
    if p.reverse,
        if any(zy(dim+1:nd)>1)
            ndxs = [ 1 : dim - 1, dim + 1 : nd, dim ];
        end
    else
        if any(zy(1:dim-1)>1)
            ndxs = [ dim, 1 : dim - 1, dim + 1 : nd ];
        end
    end
    if ~isempty( ndxs ),
        if ox,
            ndxs = permute( reshape( 1 : nx * nv, zx ), ndxs );
            x.basis_ = x.basis_( :, ndxs );
        else
            x = permute( x, ndxs );
        end
        if p.reduce, 
            ndxs = []; 
        end
    end
end
if p.reverse,
    sx = [ nv, nx ];
    xdim = 2;
    vdim = 1;
else
    sx = [ nx, nv ];
    xdim = 1;
    vdim = 2;
end

% Determine expression types
x = reshape( x, sx );
map = reshape( p.map( :, cvx_classify( x ) ), [], sx(1), sx(2) );
[ mmm, map ] = max( all( map, xdim + 1 ), [], 1 );

if ~all( mmm ),
    % Errors found
    mmm = ~mmm(:);
    sx( vdim ) = nnz( mmm );
    mmm = repmat( mmm, [ 1, nx ] );
    if ~p.reverse, mmm = mmm'; end
    x = reshape( cvx_subsref( x, mmm ), sx );
    if p.reverse, x = x'; end
    if isfield( p, 'errargs' ), varargin = varargin(p.errargs); end
    cvx_dcp_error( p.name, 'reduce', x, varargin{:} );
end

if length( mmm ) ~= 1,
    mu = sort( map(:)' ); %#ok
    mu = mu([true,diff(mu)~=0]);
else
    mu = map;
end

if length(mu) == 1,
    
    % Homogeneous input (single type)
    sx = any( mu == p.constant );
    if sx, x = cvx_constant( x ); end
    y = p.funcs{mu}( x, varargin{:} );
    if sx && ox, y = cvx( y ); end
    
else
    
    % Heterogenous input (multiple types)
    sy = size( x );
    rep = sy;
    if p.reduce, 
        sy( xdim ) = 1; 
    end
    rep( vdim ) = 1;
    if ox,
        y = cvx( sy, [] );
    else
        y = zeros( sy );
    end
    map = map(:);
    for kk = mu
        tt = map == mu;
        t2 = repmat( tt, rep );
        sy( vdim ) = nnz( tt );
        b = cvx( sz, x.basis_( :, t2 ) );
        sx = any( mu == p.constant );
        if sx, b = cvx_constant( b ); end
        b = p.funcs{kk}( b, varargin{:} );
        if sx && ox, b = cvx( b ); end
        if ~p.reduce, tt = t2; end
        if ox,
            y.basis_( 1:size(b,1), tt ) = b.basis_;
        else
            y( tt ) = b;
        end
    end
    
end

% Reverse permute/reshape (only needed if we didn't reduce)
if ~isempty( ndxs ),
    if ox,
        y.basis_(:,ndxs) = y.basis_;
    else
        y = ipermute( y, ndxs );
    end
end
y = reshape( y, zy );

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
