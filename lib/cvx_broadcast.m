function [ x, y, sz, xs, ys, sdp_mode ] = cvx_broadcast( x, y, sdp_mode, op )

global cvx___
sx = size( x ); 
sy = size( y ); 
if all( sx == 1 ),
    sz = sy;
    if all( sy == 1 ),
        xs = false;
        ys = false;
        sdp_mode = false;
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
else
    nd = max( length( sx ), length( sy ) );
    sx( end + 1 : nd ) = 1;
    sy( end + 1 : nd ) = 1;
    if any( sx ~= sy & sx ~= 1 & sy ~= 1 ),
        error( 'Matrix dimensions must agree.');
    elseif ~cvx___.broadcast,
        % Put a warning in here at some point
        error( 'Matrix dimensions must agree.');
    end
    bx = ones( 1, nd ); by = bx;
    tx = sx == 1; bx( tx ) = sy( tx );
    ty = sy == 1; by( ty ) = sx( ty );
    x = repmat( x, bx );
    y = repmat( y, by );
    sz = size( x );
end
if nargin >= 3 && sdp_mode,
    sdp_mode = false;
    if all( sz( 1 : 2 ) > 1 ),
        if sz( 1 ) ~= sz( 2 ),
            error( 'SDP constraints must be square.' );
        elseif xs && cvx_isnonzero( x ),
            error( 'SDP constraint {scalar} %s {matrix} valid only if the scalar is zero.', op );
        elseif ys && cvx_isnonzero( y ), 
            error( 'SDP constraint {matrix} %s {scalar} valid only if the scalar is zero.', op );
        elseif ~cvx_isaffine( x ) || ~cvx_isaffine( y ),
            error( 'Both sides of an SDP constraint must be affine.' );
        else
            sdp_mode = true;
        end
    end
    if ~sdp_mode,
        x = x(:);
        y = y(:);
    end
else
    x = x(:);
    y = y(:);
end



