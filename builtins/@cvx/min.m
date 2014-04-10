function z = min( x, y, dim )

%   Disciplined convex/geometric programming information:
%       MAX is convex, log-log-concave, and nondecreasing in its first
%       two arguments. Thus when used in disciplined convex programs, 
%       both arguments must be concave (or affine). In disciplined 
%       geometric programs, both arguments must be log-concave/affine.

persistent remap remap2 funcs funcs2
if nargin == 2,
    
    %
    % min( X, Y )
    %

    if isempty( remap ),
        remap = cvx_remap( ...
            { { 'real' } }, ...
            { { 'nonpositive', 'n_nonconst' }, { 'nonnegative', 'p_nonconst' } }, ...
            { { 'nonnegative', 'p_nonconst' }, { 'nonpositive', 'n_nonconst' } }, ...
            { { 'l_concave', 'positive' } }, ...
            { { 'concave' } } );
        funcs = { @min_b1, @min_b2, @min_b3, @min_b4, @min_b5 };
    end
    
    try
        z = binary_op( 'min', funcs, remap, x, y );
    catch exc
        if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
        else rethrow( exc ); end
    end

elseif nargin > 1 && ~isempty( y ),

    error( 'MIN with two matrices to compare and a working dimension is not supported.' );
        
else
    
    %
    % min( X, [], dim )
    %

    persistent params
    if isempty( params ),
        params.map     = cvx_remap( { 'real' ; 'l_concave' ; 'concave' } );
        params.funcs   = { @min_r1, @min_r2, @min_r3 };
        params.zero    = [];
        params.reduce  = true;
        params.reverse = false;
        params.dimarg  = 2;
        params.name    = 'min';
    end
    
    if nargin > 1 && ~isempty( y ),
        error( 'MIN with two matrices to compare and a working dimension is not supported. ');
    end
    if nargin < 3, 
        dim = []; 
    end
    
    try
        y = reduce_op( params, x, dim );
    catch exc
        if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
        else rethrow( exc ); end
    end
   
end 

function z = min_b1( x, y )
% constant
z = cvx( min( cvx_constant( x ), cvx_constant( y ) ) );

function z = min_b2( x, y ) %#ok
% min( nonnegative, nonpositive ) (pass-through x)
z = x;

function z = min_b3( x, y ) %#ok
% min( nonpositive, nonnegative ) (pass-through y)
z = y;

function z = min_b4( x, y )
z = [];
sz = [max(numel(x),numel(y)),1]; %#ok
% geometric
cvx_begin gp
    hypograph variable z( sz )
    x >= z; %#ok
    y >= z; %#ok
cvx_end

function z = min_b5( x, y )
z = [];
sz = [max(numel(x),numel(y)),1]; %#ok
% linear
xsg = cvx_sign( x );
ysg = cvx_sign( y );
xzr = ~cvx_isnonzero( x, true );
yzr = ~cvx_isnonzero( y, true );
cvx_begin
    hypograph variable z( sz )
    x >= z; %#ok
    y >= z; %#ok
    cvx_setnpos( fastref( z, (xsg<0)|xzr|(ysg<0)|yzr ) );
    cvx_setnneg( fastref( z, ((xsg>0)|xzr)&((ysg>0)|yzr) ) );
cvx_end

function x =  min_r1( x )
x = min( cvx_constant( x ), [], 1 );

function x = min_r2( x )
nx = x.size_(1);
if nx > 1,
    nv = x.size_(2); %#ok
    cvx_begin gp
        hypograph variable z( 1, nv )
        x >= repmat( z, nx, 1 ); %#ok
    cvx_end
    x = cvx_optval;
end

function x = min_r3( x )
nx = x.size_(1);
if nx > 1,
    nv = x.size_(2); %#ok
    xsg = cvx_sign( xt );
    xzr = ~cvx_isnonzero( xt, true );
    cvx_begin
        hypograph variable zt( 1, nv )
        xt >= repmat( zt, nx, 1 ); %#ok
        cvx_setnpos( fastref( zt, any((xsg<0)|xzr,1) ) );
        cvx_setnneg( fastref( zt, all((xsg>0)|xzr,1) ) );
    cvx_end
    x = cvx_optval;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
