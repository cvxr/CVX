function [ p, pstr ] = cvx_validate( p, id, clean )
global cvx___
np = length( cvx___.problems );
if p <= np,
    pstr = cvx___.problems( p );
    if pstr.id == id, 
        if p > np && ( nargin < 3 || clean ),
            cvx_pop( p + 1, true );
        end
        return
    end
end
cvx_pop( 0, true );
cvx_throw( 'Internal CVX data corruption. Please CLEAR ALL and rebuild your model.' );

