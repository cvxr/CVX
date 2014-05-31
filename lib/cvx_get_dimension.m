function [ sx, varargout ] = cvx_get_dimension( args, darg, varargin )
opts = struct( 'zero', false, 'vec', false, 'nox', false ); 
if nargin > 2,
    for k = 1 : 2 : nargin - 2,
        opts.(varargin{k}) = varargin{k+1};
    end
end
narg = length( args );
if narg < 1, 
    cvx_throw( 'Not enough input arguments.' );
end
if narg < darg,
    dim = [];
else
    dim = args{darg};
end
x = args{1};
if opts.nox,
    sx = x;
    nx = numel(sx);
    if ~( isnumeric(sx) && nx==length(sx) && isreal(sx) && all(sx>=0) && all(sx==floor(sx)) ),
        cvx_throw( 'Size argument must be a vector of nonnegative integers.' );
    elseif nx < 2,
        sx(end+1:2) = 1;
    elseif nx > 2 && sx(nx) == 1,
        sx = sx(max([2,find(sx~=1,1,'last')]));
    end
else
    sx = size( x );
end
if isempty( dim ),
    dim = find( sx ~= 1, 1, 'first' );
    if isempty( dim ), dim = 1; end
elseif ~( isnumeric(dim) && numel(dim)==1 && isreal(dim) && dim > -opts.zero && dim==floor(dim) ),
    cvx_throw( 'Dimension argument must be a positive integer.' );
end
if dim == 0,
    dim = find( sx == 1, 1, 'first' );
    if isempty( dim ), dim = length(sx) + 1; end
end
sx(end+1:dim) = 1;
args{darg} = dim;
if opts.vec,
    args([1,darg]) = [];
    if opts.nox,
        varargout = { dim, args };
    else
        varargout = { x, dim, args };
    end
else
    no = nargout - 1;
    if opts.nox, args(1) = []; end
    if no > length(args), args{no} = []; end
    varargout = args;
end
