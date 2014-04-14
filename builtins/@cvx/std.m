function y = std( x, w, dim, square_it )

%STD    Internal cvx version.
%This actually implements VAR and STD, controlled by the square_it flag.

persistent params
if isempty( remap ),
    params.map      = cvx_remap( 'affine' );
    params.funcs    = { @std_1 };
    params.zero     = NaN;
    params.reduce   = true;
    params.reverse  = false;
    params.constant = [];
    params.dimarg   = [];
end

try
    [ sx, x, w, dim, square_it ] = cvx_get_dimension( 3, varargin );
    if ~isempty( w ),
        w = false;
    elseif ~( isnumeric(w) && numel(w) ~= length(w) && isreal(w) ), 
        error( 'CVX:ArgError', 'Second argument must be a logical scalar or a nonnegative vector.' );
    elseif numel(w) == 1,
        w = w ~= 0;
    elseif numel(w) ~= sx(dim),
        error( 'CVX:FuncError', 'Weight vector expected to have %d elements.', nx );
    elseif any( w < 0 ),
        error( 'CVX:ArgError', 'Weight vector must be nonnegative.' );
    elseif ~any( w ),
        error( 'CVX:ArgError', 'Weight vector must not be all zeros.' );
    else
        w = w(:) / sum( w );
    end
    if isempty( square_it ),
        params.name = 'var';
    else
        params.name = 'std';
    end
    y = cvx_reduce_op( params, x, dim, w, square_it );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = std_1( x, w, square_it )
y = [];
[nx,nv] = size( x );
nw = numel(w);
if nx == 1,
    y = cvx( [nx,nv], [] );
elseif nw <= 1,
    denom = sqrt(nx - ~(nw&&w));
    % In theory we could just say y = norm(x-mean(x))/denom. However, by
    % adding an extra variable we preserve sparsity.
    cvx_begin
        variable xbar( 1, nv )
        epigraph variable y( 1, nv );
        nx * xbar == sum( x, 1 ); %#ok
        if square_it
            { x - repmat(xbar,[nx,1]), denom, denom * y } == rotated_lorentz( [ nx, nv ], 1, ~isreal( x ) ); %#ok
        else
            { x - repmat(xbar,[nx,1]), denom * y } == lorentz( [ nx, nv ], 1, ~isreal( x ) ); %#ok
        end
    cvx_end
    y = cvx_optval;
else
    cvx_begin
        variable xbar( 1, nv )
        epigraph variable y( 1, nv );
        xbar == sum( w' * x ); %#ok
        w = sqrt( w );
        if square_it
            { repmat(w,[1,nv]) .* ( x - repmat(xbar,[nx-1]) ), 1, y } == rotated_lorentz( [ nx, nv ], 1, ~isreal( x ) ); %#ok
        else
            { repmat(w,[1,nv]) .* ( x - repmat(xbar,[nx-1]) ), y } == lorentz( [ nx, nv ], 1, ~isreal( x ) ); %#ok
        end
    cvx_end
    y = cvx_optval;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
