function [ sx, x, dim, varargout ] = cvx_get_dimension( darg, args, vec ) %#ok

narg = length( args );
oarg = nargout - 1;
if narg < 1, 
    error( 'CVX:ArgError', 'Not enough input arguments.' );
end
zero_ok = darg < 0;
if zero_ok,
    darg = -darg;
end
if narg < darg,
    dim = [];
else
    dim = args{darg};
end
x = args{1};
sx = size( x );
if isempty( dim ),
    dim = find( sx ~= 1, 1, 'first' );
    if isempty( dim ), dim = 1; end
elseif ~( isnumeric(dim) && numel(dim)==1 && isreal(dim) && dim > -zero_ok && dim==floor(dim) ),
    error( 'CVX:ArgError', 'Dimension argument must be a positive integer.' );
end
if dim == 0,
    dim = find( sx == 1, 1, 'first' );
    if isempty( dim ), dim = length(sx)+1; end
end
sx(end+1:dim) = 1;
args{darg} = dim;
if nargin < 3,
    args{oarg+1} = [];
    varargout = args(3:end);
else
    args([1,darg]) = [];
    varargout{1} = args;
end
