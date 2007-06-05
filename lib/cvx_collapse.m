function x = cvx_collapse( x, keeptemp, tocell )
if nargin < 2, keeptemp = false; end
if nargin < 3, tocell = false; end

while true,
    sx = size( x );
    nx = prod( sx );
    switch class( x ),
        case 'cell',
            if nx == 1,
                x = x{1};
                continue;
            end
            x = reshape( x, 1, nx );
        case 'struct',
            if keeptemp,
                fx = fieldnames( x );
            else
                [ fx, ndxs ] = cvx_fieldnames( x );
            end
            nfx = length( fx );
            if nfx == 1 & nx == 1,
                x = subsref( x, struct( 'type', '.', 'subs', fx{1} ) );
                continue;
            end
            if tocell,
                if nfx == 1,
                    sx = [ 1, sx ];
                else
                    sx = [ 1, nfx, sx ];
                end
                x = struct2cell( x );
                if ~keeptemp,
                    x = x( ndxs, : );
                end
                x = reshape( x, sx );
            end
        otherwise,
            cx = false;
    end
    break;
end
