function y = cvx_reduce_op( p, varargin )

% Get dimension and main argument
[ zx, x, dim, args ] = cvx_get_dimension( varargin, p.dimarg, 'vec', true );
ox = isa( x, 'cvx' );
nd = length( zx );

% Quick exit for empty results
nx = zx( dim );
zy = zx;
if p.reduce && ( nx || ~isempty( p.zero ) ),
    zy( dim ) = 1; 
end
if ~all( zy ),
    if ~any( zy ), zy = [ 1, 1 ]; end
    y = cvx( repmat( p.zero ), zy );
    return
end
nv = prod( zx ) / nx;

% Permute if needed, and reshape to canonical size
perm = [];
if nx > 1 && nv > 1,
    if p.reverse,
        if any(zy(dim+1:nd)>1)
            perm = [ 1 : dim - 1, dim + 1 : nd, dim ];
        end
    else
        if any(zy(1:dim-1)>1)
            perm = [ dim, 1 : dim - 1, dim + 1 : nd ];
        end
    end
    if ~isempty( perm ),
        if ox,
            x    = cvx_basis( x );
            ndxs = permute( reshape( 1 : nx * nv, zx ), perm );
            x    = cvx( zx(perm), x(:,ndxs) );
        else
            x = permute( x, perm );
        end
        if p.reduce, 
            perm = [];
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
if isempty( p.map ),
    mu = 1;
else
    map = cvx_classify( x );
    map = p.map( :, map );
    map = reshape( map, [], sx(1), sx(2) );
    map = all( map, xdim + 1 );
    [ mmm, map ] = max( map, [], 1 );
    if ~all( mmm ),
        mu = 0;
    elseif length( mmm ) ~= 1,
        mu = sort( map(:)' ); %#ok
        mu = mu([true,diff(mu)~=0]);
    else
        mu = map;
    end
end

if length(mu) == 1,
    
    % Homogeneous input (single type)
    if mu ~= 0,
        sx = any( mu == p.constant );
        if sx, x = cvx_constant( x ); end
        y = p.funcs{mu}( x, args{:} );
        if sx && ox, y = cvx( y ); end
    end
    
    % Post-op check, if necessary
    if isempty( p.map ),
        mu = cvx_isvalid( y );
        if ~mu, 
            mmm = cvx_isvalid( y, true ); 
        end
    end
        
    % Errors found
    if mu == 0,
        mmm = ~mmm(:);
        sx( vdim ) = nnz( mmm );
        mmm = repmat( mmm, [ 1, nx ] );
        if ~p.reverse, mmm = mmm'; end
        x = reshape( cvx_subsref( x, mmm ), sx );
        if p.reverse, x = x'; end
        if isfield( p, 'errargs' ), args = args(p.errargs); end
        cvx_dcp_error( p.name, 'reduce', x, args{:} );
    end

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
    if ~p.reverse,
        map = map';
    end
    for kk = mu
        tt = map == kk;
        sz = [ nx, nnz( tt ) ];
        t2 = repmat( tt, rep );
        b = reshape( cvx_subsref( x, t2 ), sz );
        sx = any( kk == p.constant );
        if sx, b = cvx_constant( b ); end
        b = p.funcs{kk}( b, args{:} );
        if sx && ox, b = cvx( b ); end
        if ~p.reduce, tt = t2; end
        if ox,
            y = cvx_subsasgn( y, tt, b );
        else
            y( tt ) = b;
        end
    end
    
end

% Reverse permute/reshape (only needed if we didn't reduce)
if isempty( perm ),
    y = reshape( y, zy );
elseif ox,
    y = cvx_basis( y );
    y( :, ndxs ) = y;
    y = cvx( zy, y );
else
    y = reshape( ipermute( y, perm ), zy );
end

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
