function [ p, pstr ] = verify( prob, clean )
global cvx___
p = prob.index_;
np = length( cvx___.problems );
if p == np,
    pstr = cvx___.problems( p );
    if cvx_id( pstr.self ) == cvx_id( prob ), return; end
elseif p < np,
    pstr = cvx___.problems( p );
    if cvx_id( pstr.self ) == cvx_id( prob ),
        if nargin < 2 || clean, cvx_pop( p + 1, true ); end
        return; 
    end
end
cvx_pop( 0 );
error( 'CVX:InternalError', 'Internal CVX data corruption. Please CLEAR ALL and rebuild your model.' );
