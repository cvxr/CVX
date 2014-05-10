function ndim = cvx_newvar( nx, type )

if nx == 0,
    ndim = zeros(0,1);
    return
end
global cvx___
persistent aff_code
if isempty( aff_code ),
    aff_code = int8(find(cvx_remap('affine_')));
end
if nargin < 2,
    type = aff_code;
end
nO = length( cvx___.classes );
ndim = nO + 1 : nO + nx;
cvx___.classes ( ndim, 1 ) = type(:);
