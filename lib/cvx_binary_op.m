function [ z, sdp_mode ] = cvx_binary_op( p, x, y, varargin )

global cvx___

sx = size( x );
sy = size( y );
if ~isempty( p.map ),
    vx = cvx_classify( x );
    vy = cvx_classify( y );
end

% Determine if sizes are compatible or broadcastable
if all( sx == 1 ),
    sz = sy;
    if all( sy == 1 ),
        xs = false;
        ys = false;
    else
        xs = true;
        ys = false;
    end
elseif all( sy == 1 ),
    sz = sx;
    ys = true;
    xs = false;
elseif isequal( sx, sy );
    sz = sx;
    xs = false;
    ys = false;
elseif ~cvx___.broadcast,
    if length( sx ) < 2 && length( sy ) < 2,
        cvx_throw( 'Matrix dimensions must agree.' );
    else
        cvx_throw( 'Array dimensions must match for binary array op.' );
    end
else
    nd = max( length( sx ), length( sy ) );
    sx( end + 1 : nd ) = 1;
    sy( end + 1 : nd ) = 1;
    if any( sx ~= sy & sx ~= 1 & sy ~= 1 ),
        cvx_throw( 'Non-singleton dimensions of the two input arrays must match each other.' );
    end
    bx = ones( 1, nd ); by = bx;
    tx = sx == 1; bx( tx ) = sy( tx );
    ty = sy == 1; by( ty ) = sx( ty );
    x = repmat( x, bx );
    y = repmat( y, by );
    vx = repmat( reshape( vx, sx ), bx );
    vy = repmat( reshape( vy, sy ), by );
    sz = size( x );
    xs = false;
    ys = false;
end

% Special size checks for SDP mode
sdp_mode = false;
if isfield( p, 'sdp' ) && p.sdp && ~( xs && ys ),
    if all( sz( 1 : 2 ) > 1 ),
        persistent sdpmap %#ok
        if isempty( sdpmap ),
            sdpmap = cvx_remap( 'affine' );
        end
        if sz( 1 ) ~= sz( 2 ),
            cvx_throw( 'SDP constraints must be square.' );
        elseif xs && cvx_isnonzero( x ),
            cvx_throw( 'SDP constraint {scalar} %s {matrix} valid only if the scalar is zero.', op );
        elseif ys && cvx_isnonzero( y ), 
            cvx_throw( 'SDP constraint {matrix} %s {scalar} valid only if the scalar is zero.', op );
        elseif ~all(sdpmap(vx(:))) || ~all(sdpmap(vy(:))),
            cvx_throw( 'SDP constraints must be affine.' );
        else
            sdp_mode = true;
        end
    end
end

% Empty arguments
if ~all( sz ),
    if isa( x, 'cvx' ) || isa( y, 'cvx' ),
        z = cvx( sz, [] );
    else
        z = zeros( sz );
    end
    return
end

% Normal expression maps
if sdp_mode || isempty( p.map ),
    vu = 1;
else
    vz = int32(vx) + size(p.map,1) * int32(vy-1);
    vz = p.map( vz );
    vu = sort( vz(:)' ); %#ok
    vu = vu( [ true, diff(vu) ~= 0 ] );
end

if vu(1) ~= 0,
    if length( vu ) == 1,
        % Homogenous input (single compute mode)
        if any( vu == p.constant ),
            x = cvx_constant( x );
            y = cvx_constant( y );
            oz = isa( x, 'cvx' ) || isa( y, 'cvx' );
        else
            oz = false;
        end  
        z = p.funcs{vu(1)}( vec(x), vec(y), varargin{:} );
        % Post-op check, if necessary
        if isempty( p.map ),
            vu = cvx_isvalid( z );
            if vu == 0, 
                vz = cvx_isvalid( z, true ); 
            end
        end
        % Output CVX object even for constant data, if requested
        z = reshape( z, sz );
        if oz, 
            z = cvx( z ); 
        end
    else
        % Heterogeneous input (multiple compute modes)
        if isa( x, 'cvx' ) || isa( y, 'cvx' ),
            z = cvx( sz, [] );
        else
            z = zeros( sz );
        end
        if xs, xt = x; end
        if ys, yt = y; end
        for vk = vu,
            tt = vz == vk;
            if ~xs, xt = cvx_fastref( x, tt ); end
            if ~ys, yt = cvx_fastref( y, tt ); end
            if any( vk == p.constant ),
                xt = cvx_constant( x );
                yt = cvx_constant( y );
            end
            zt = p.funcs{vk}( xt, yt, varargin{:} );
            z = cvx_fastasgn( z, tt, zt );
        end
    end
end

% Errors found
if vu(1) == 0,
    if ~xs, x = cvx_fastref( x, vz == 0 ); end
    if ~ys, y = cvx_fastref( y, vz == 0 ); end
    cvx_dcp_error( p.name, 'binary', x, y );
end



    
