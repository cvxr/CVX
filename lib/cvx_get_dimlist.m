function varargout = cvx_get_dimlist( args, varargin )
opts = struct( 'default', [], 'nd', true );
if nargin > 2,
    for k = 1 : 2 : nargin - 2,
        opts.(varargin{k}) = varargin{k+1};
    end
end
narg = length( args );
if narg < 1,
    if isempty( opts.default ),
        error( 'CVX:ArgError', 'Not enough input arguments.' );
    end
    sx = [];
else
    sx = args{1};
end
if isempty( sx ),
    sx = opts.default;
    if isempty( sx ),
        error( 'CVX:ArgError', 'Size argument must be a vector of nonnegative integers.' );
    end
elseif ~( isnumeric(sx) && numel(sx)==length(sx) && isreal(sx) && all(sx>=0) && all(sx==floor(sx)) ),
    error( 'CVX:ArgError', 'Size argument must be a vector of nonnegative integers.' );
end
if length( sx ) > 2,
    n = max( [2,find( sx > 1, 1, 'last' )]);
    sx(n+1:end) = [];
    if ~opts.nd && length( sx ) > 2,
        error( 'CVX:ArgError', 'N-D arrays are not supported.' );
    end
else
    sx(end+1:2) = 1;
end
args{1} = sx;
if nargout > length(args), args{nargout} = []; end
varargout = args;
