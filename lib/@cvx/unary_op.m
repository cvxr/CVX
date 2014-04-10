function y = unary_op( fname, funcs, remap, x, varargin )
sx = size( x );
if ~all( sx ),
    y = x;
    return
end
vx = remap( cvx_classify( x ) );
vu = sort( vx(:)' ); %#ok
if vu( 1 ) == 0,
    if ~isempty(varargin) && isnumeric(varargin{1}) && numel(varargin{1}) == 1,
        cvx_dcp_error( x, varargin{1}, vx == 0, fname );
    else
        cvx_dcp_error( x, vx == 0, fname );
    end
end
vu = vu( [ true, diff(vu) ~= 0 ] );
if length( vu ) == 1,
    x.size_ = [prod(sx),1];
    y = funcs{vu(1)}( x, varargin{:} );
    y.size_ = sx;
else
    y = cvx( sx, [] ); 
    for vk = vu,
        tt = vx == vk;
        b = x.basis_( :, tt );
        b = cvx( [size(b,2),1], b );
        b = funcs{vk}( b, varargin{:} );
        b = b.basis_;
        y.basis_(1:size(b,1),tt) = b;
    end
end
