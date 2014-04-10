function z = binary_op( fname, funcs, remap, x, y, varargin )

[ x, y, sz, xs, ys ] = cvx_broadcast( x, y );
if ~all( sz ),
    z = cvx( sz, [] );
    return
end

cx = isa( x, 'cvx' );
cy = isa( y, 'cvx' );
vx = cvx_classify( x );
vy = cvx_classify( y );
vz = remap( vx + size( remap, 1 ) * ( vy - 1 ) );
vu = sort( vz(:)' ); %#ok
if vu( 1 ) == 0,
    cvx_dcp_error( x, y, vz == 0, fname );
end
vu = vu( [ true, diff(vu) ~= 0 ] );
if length( vu ) == 1,
    z = funcs{vu(1)}( x, y, varargin{:} );
    z.size_ = sz;
else
    z = cvx( sz, [] );
    if xs, xt = x; end
    if ys, yt = y; end
    for vk = vu,
        tt = vz == vk;
        if ~xs,
            if cx,
                b = x.basis_( :, tt );
                xt = cvx( [size(b,2),1], b );
            else
                xt = x(tt);
            end
        end
        if ~ys,
            if cy,
                b = y.basis_( :, tt );
                yt = cvx( [size(b,2),1], b );
            else
                yt = y(tt);
            end
        end
        zt = funcs{vk}( xt, yt, varargin{:} );
        b = zt.basis_;
        z.basis_(1:size(b,1),tt) = b;
    end
end
