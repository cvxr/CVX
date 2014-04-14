function y = cvx_unary_op( p, x, varargin )
sx = size( x );
ox = isa( x, 'cvx' );
if ~all( sx ),
    y = x;
    return
end
vx = p.map( cvx_classify( x ) );
vu = sort( vx(:)' ); %#ok
if vu( 1 ) == 0,
    x = cvx_subsref( x, vx == 0 );
    cvx_dcp_error( p.name, 'unary', x, varargin{:} );
end
vu = vu( [ true, diff(vu) ~= 0 ] );
if length( vu ) == 1,
    xs = any( vu == p.constant );
    if xs, x = cvx_constant( x ); end
    y = p.funcs{vu(1)}( x(:), varargin{:} );
    y = reshape( y, sx );
    if ox, y = cvx( y ); end
elseif all( p.constant(vu) ),
    if ox, y = cvx( sx, [] );
    else y = zeros( sx ); end
    for vk = vu,
        tt = vx == vk;
        xs = any( vk == p.constant );
        b = vec( cvx_subsref( x, tt ) );
        if xs, b = cvx_constant( b ); end
        b = p.funcs{vk}( b, varargin{:} );
        if ox, b = cvx( b ); end
        y = cvx_subsasgn( y, tt, b );
    end
end
