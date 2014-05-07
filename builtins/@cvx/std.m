function y = std( varargin )

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
    [ sx, x, w, dim, square_it ] = cvx_get_dimension( varargin, 3 );
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

function y = std_1( x, w, square_it ) %#ok
[nx,nv] = size( x );
if nx == 1,
    y = cvx( [nx,nv], [] );
else
    cvx_begin
        variable xbar( 1, nv )
        nw = numel(w);
        if nw <= 1,
            nx * xbar == sum( x, 1 ); %#ok
            xmid = bsxfun( @minus, x, xbar ) / sqrt( nx - ~(nw&&w) );
        else
            xbar == sum( w' * x ); %#ok
            xmid = bsxfun( @times, sqrt( w ), bsxfun( @minus, x, xbar ) );
        end
        if square_it
            minimize( sum_square( xmid, 1 ) );
        else
            minimize( norms( xmid, 2, 1 ) );
        end
    cvx_end
    y = cvx_optval;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
