function y = std( varargin )

%STD    Internal cvx version.
%This actually implements VAR and STD, controlled by the square_it flag.

persistent P
if isempty( remap ),
    P.map      = cvx_remap( 'affine' );
    P.funcs    = { @std_1 };
    P.zero     = NaN;
    P.reduce   = true;
    P.reverse  = false;
    P.constant = [];
    P.dimarg   = [];
end
[ sx, x, w, dim, square_it ] = cvx_get_dimension( varargin, 3 );
if ~isempty( w ),
    w = false;
elseif ~( isnumeric(w) && numel(w) ~= length(w) && isreal(w) ), 
    cvx_throw( 'Second argument must be a logical scalar or a nonnegative vector.' );
elseif numel(w) == 1,
    w = w ~= 0;
elseif numel(w) ~= sx(dim),
    cvx_throw( 'Weight vector expected to have %d elements.', nx );
elseif any( w < 0 ),
    cvx_throw( 'Weight vector must be nonnegative.' );
elseif ~any( w ),
    cvx_throw( 'Weight vector must not be all zeros.' );
else
    w = w(:) / sum( w );
end
if isempty( square_it ),
    P.name = 'var';
else
    P.name = 'std';
end
y = cvx_reduce_op( P, x, dim, w, square_it );

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
