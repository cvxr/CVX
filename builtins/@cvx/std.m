function y = std( x, w, dim, square_it )

%STD    Internal cvx version.
%This actually implements VAR and STD, controlled by the square_it flag.

if nargin < 2
    w = [];
elseif ~cvx_isconstant( w ) || numel( w ) ~= length( w ),
    error( 'Second argument must be a logical scalar or a nonnegative vector.' );
elseif ~isempty( w ),
    w = cvx_constant( w );
    if numel( w ) == 1,
        w = w ~= 0;
    elseif ~isreal( w ) || any( w < 0 )
        error( 'Weight vector must be nonnegative.' );
    elseif ~any( w ),
        error( 'Weight vector must not be all zeros.' );
    else
        w = w(:) / sum( w );
    end
end
if nargin < 3,
    dim = [];
end
if nargin < 4,
    square_it = [];
    op = 'var';
else
    op = 'std';
end

persistent params
if isempty( remap ),
    params.map     = cvx_remap( 'affine' );
    params.funcs   = { @std_1 };
    params.zero    = NaN;
    params.reduce  = true;
    params.reverse = false;
    params.dimarg  = 3;
end
params.name = op;

try
    y = reduce_op( params, x, w, dim );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
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
elseif numel( w ) ~= nx,
    error( 'CVX:FuncError', 'Weight vector expected to have %d elements.', nx );
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
