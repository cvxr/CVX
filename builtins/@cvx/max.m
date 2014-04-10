function z = max( x, y, dim )

%   Disciplined convex/geometric programming information:
%       MAX is convex, log-log-convex, and nondecreasing in its first 
%       two arguments. Thus when used in disciplined convex programs, 
%       both arguments must be convex (or affine). In disciplined 
%       geometric programs, both arguments must be log-convex/affine.

persistent remap remap2 funcs funcs2
if nargin == 2,
    
    %
    % max( X, Y )
    %

    if isempty( remap ),
        remap = cvx_remap( ...
            { { 'real' } }, ...
            { { 'nonnegative', 'p_nonconst' }, { 'nonpositive', 'n_nonconst' } }, ...
            { { 'nonpositive', 'n_nonconst' }, { 'nonnegative', 'p_nonconst' } }, ...
            { { 'l_convex', 'positive' } }, ...
            { { 'convex' } } );
        funcs = { @max_b1, @max_b2, @max_b3, @max_b4, @max_b5 };
    end
    
    try
        z = binary_op( 'max', funcs, remap, x, y );
    catch exc
        if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
        else rethrow( exc ); end
    end

elseif nargin > 1 && ~isempty( y ),

    error( 'MAX with two matrices to compare and a working dimension is not supported.' );
        
else
    
    %
    % max( X, [], dim )
    %

    persistent params
    if isempty( params ),
        params.map     = cvx_remap( { 'real' ; 'l_convex' ; 'convex' } );
        params.funcs   = { @max_r1, @max_r2, @max_r3 };
        params.zero    = [];
        params.reduce  = true;
        params.reverse = false;
        params.dimarg  = 2;
        params.name    = 'max';
    end

    if nargin > 1 && ~isempty( y ),
        error( 'MAX with two matrices to compare and a working dimension is not supported. ');
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

function z = max_b1( x, y )
% constant
z = cvx( max( cvx_constant( x ), cvx_constant( y ) ) );

function z = max_b2( x, y )
% max( positive, negative ) (pass through x)
z = x;

function z = max_b3( x, y )
% max( negative, positive ) (pass through y)
z = y;

function z = max_b4( x, y )
% geometric
z = [];
sz = [ max(numel(x),numel(y)), 1 ]; %#ok
cvx_begin gp
    epigraph variable z( sz );
    x <= z; %#ok
    y <= z; %#ok
cvx_end

function z = max_b5( x, y )
% linear
z = [];
sz = [ max(numel(x),numel(y)), 1 ]; %#ok
xsg = cvx_sign( x );
ysg = cvx_sign( y );
xzr = ~cvx_isnonzero( x, true );
yzr = ~cvx_isnonzero( y, true );
cvx_begin
    epigraph variable z( sz );
    x <= z; %#ok
    y <= z; %#ok
    cvx_setnneg( fastref( z, (xsg>0)|(ysg>0)|xzr|yzr ) );
    cvx_setnpos( fastref( z, ((xsg<0)|xzr)&((ysg<0)|yzr) ) );
cvx_end

function x = max_r1( x )
if x.size_ > 1,
    x = max( cvx_constant( x ), [], 1 );
end

function x = max_r2( x )
nx = x.size_(1);
if nx > 1,
    nv = x.size_(2); %#ok
    cvx_begin gp
        epigraph variable z( 1, nv )
        x <= repmat( z, nx, 1 ); %#ok
    cvx_end
    x = cvx_optval;
end

function x = max_r3( x )
nx = x.size_(1);
if nx > 1,
    nv = x.size_(2); %#ok
    xsg = cvx_sign( x );
    xzr = ~cvx_isnonzero( x, true );
    cvx_begin
        epigraph variable z( 1, nv )
        x <= repmat( z, nx, 1 ); %#ok
        cvx_setnneg( fastref( z, any((xsg>0)|xzr,1) ) );
        cvx_setnpos( fastref( z, all((xsg<0)|xzr,1) ) );
    cvx_end
    x = cvx_optval;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
