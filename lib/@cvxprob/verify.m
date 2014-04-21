function [ p, pstr ] = verify( prob, clean )
global cvx___
p = prob.index_;
try
    np = length( cvx___.problems );
catch
    np = 0;
end
if p <= np,
    pstr = cvx___.problems( p );
end
if p > np || pstr.self ~= prob,
    cvx_pop( 0 );
    error( 'CVX:InternalError', 'Internal CVX data corruption. Please CLEAR ALL and rebuild your model.' );
elseif np < p && ( nargin < 2 || clean )
    cvx_pop( p + 1, true );
end
