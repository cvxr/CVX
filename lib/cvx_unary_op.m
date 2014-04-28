function y = cvx_unary_op( p, x, varargin )

sx = size( x );
if ~all( sx ),
    y = x;
    return
end

% Expression maps
if isempty( p.map ),
    vu = -1;
else
    vx = p.map( cvx_classify( x ) );
    vu = sort( vx(:)' ); %#ok
    vu = vu( [ true, diff(vu) ~= 0 ] );
end

if vu( 1 ) ~= 0,

    ox = isa( x, 'cvx' );
    if length( vu ) == 1,

        % Homogeneous input (single compute mode)
        if ox && vu == 1, x = cvx_constant( x ); end
        y = p.funcs{vu(1)}( vec(x), varargin{:} );

        % Post-op check, if necessary
        if isempty( p.map ),
            vu = cvx_isvalid( y );
            if vu == 0,
                vx = cvx_isvalid( y, true ); 
            end
        end

        % Output CVX object even for constant data, if requested
        y = reshape( y, sx );
        if ox, y = cvx( y ); end

    else

        % Heterogeneous input (multiple compute modes)
        if ox, y = cvx( sx, [] );
        else y = zeros( sx ); end
        for vk = vu,
            tt = vx == vk;
            xt = cvx_fastref( x, tt );
            if ox && vk == 1, xt = cvx_constant( xt ); end
            yt = p.funcs{vk}( xt, varargin{:} );
            y = cvx_fastasgn( y, tt, yt );
        end

    end

end

if vu(1) == 0,
    x = cvx_fastref( x, ~vx );
    if isfield( p, 'name' ) && ~isempty( p.name ),
        name = p.name;
    else
        dd = dbstack;
        name = dd(2).name;
    end
    cvx_dcp_error( name, 'unary', x, varargin{:} );
end

