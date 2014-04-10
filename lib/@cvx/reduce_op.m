function y = reduce_op( p, x, varargin )

if nargin < 2,
    error( 'Not enough arguments.' );
end

% Determine reducing dimension
zx = size( x );
if nargin < p.dimarg + 1,
    dim = [];
else
    dim = varargin{p.dimarg-1};
    varargin(p.dimarg-1) = [];
end
if isempty( dim ),
    dim = find( zx ~= 1, 1, 'first' );
    if isempty( dim ), dim = 1; end
elseif ~( isnumeric( dim ) && numel( dim ) == 1 && isreal( dim ) && dim > 0 && dim ~= floor(dim) && ~isinf(dim) ),
    error( 'Dimension argument must be a positive integer scalar within indexing range.' );
end
nd = length( zx );
if nd < dim,
    zx( end + 1 : nd ) = 1;
    dim = nd;
end

% Quick exit for empty results
nx = zx( dim );
zy = zx;
if p.reduce && ( nx || ~isempty( p.zero ) ),
    zy( dim ) = 1; 
end
if ~all( zy ),
    y = cvx( repmat( p.zero ), zy )
    return
end
nv = prod( zx ) / nx;

% Permute if needed, and reshape to canonical size
ndxs = [];
if nx > 1,
    if ~p.reverse && any(zy(1:dim-1)>1),
        ndxs = [ dim, 1 : dim - 1, dim + 1 : nd ];
    elseif any(zy(dim+1:nd)>1),
        ndxs = [ 1 : dim - 1, dim + 1 : nd, dim ];
    end
    if isempty( ndxs ),
        ndxs = permute( reshape( 1 : nx * nv, zx ), perm );
        x.basis_ = x.basis_(:,ndxs);
        if p.reduce, ndxs = []; end
    end
end
if p.reverse,
    x.size_ = [ nv, nx ];
else
    x.size_ = [ nx, nv ];
end

map = reshape( p.map( :, cvx_classify( x ) ), [], x.size_ );
[ mmm, map ] = max( all( map, 3 - ~p.reverse ), [], 1 );
if ~all( mmm ),
    mmm = ~mmm(:);
    x.size_( 2 - ~p.reverse ) = nnz( mmm );
    mmm = repmat( mmm, [ nx, 1 ] );
    if p.reverse, mmm = mmm'; end
    x.basis_ = x.basis_( :, mmm );
    if p.reverse, x = x'; end
    cvx_dcp_error( x, p.name );
end
map = map(:)';
mu  = sort( map );
mu  = mu([true,diff(mu)~=0]);
if length(mu) == 1,
    y = p.funcs{mu}( x, varargin{:} );
else
    if p.reduce,
        y = cvx( nv, [] );
    else
        y = cvx( x.size_, [] );
    end
    for kk = mu
        tt = map == mu;
        if p.reverse,
            t2 = repmat( tt', [1,nx] );
            sz = [ nnz(tt), nx ];
        else
            t2 = repmat( tt, [nx,1] );
            sz = [ nx, nnz(tt) ];
        end
        b = cvx( sz, x.basis_( :, t2 ) );
        b = p.funcs{kk}( b, varargin{:} );
        b = b.basis_;
        if ~p.reduce, tt = t2; end
        y.basis_(1:size(b,1),tt) = b;
    end
end

% Reverse permute/reshape (only needed if we didn't reduce)
if ~isempty( ndxs ),
    y.basis_(:,ndxs) = y.basis_;
end
y.size_ = zy;

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
